function results = acoustic_indices(h, fs, bandsPerOctave, freqLims, tEarly, alignThreshold_dB)
%ACOUSTIC_INDICES  Room-acoustic indices from a mono room impulse response.
%
%   RESULTS = ACOUSTIC_INDICES(H, FS, BANDSPEROCTAVE, FREQLIMS, TEARLY,
%       ALIGNTHRESHOLD_DB)
%
%   Computes a standard set of room-acoustic parameters from a single-channel
%   room impulse response H sampled at FS Hz.
%
%   PROCESSING CHAIN
%     1. Onset alignment: the RIR is trimmed so that t = 0 corresponds
%        approximately to the direct-sound arrival, detected as the first
%        sample whose magnitude exceeds ALIGNTHRESHOLD_DB below the global
%        peak.  A small pre-offset (5 samples) is applied to capture the
%        leading edge of the direct sound.
%     2. Mid-frequency T60 crop: a quick T20 Schroeder estimate is made
%        from a single 500–1000 Hz bandpass of the aligned RIR.  The signal
%        is then hard-trimmed to that T60 length before further processing.
%        Samples beyond one T60 carry negligible energy and only inflate
%        computation.  If the estimate fails the full signal is kept.
%     3. Band-wise reverberation time using the noise-compensated Schroeder
%        integral.  The algorithm follows Lundeby et al. (1995) and the
%        ISO 3382-1:2009 Chu correction:
%          • Iteratively estimate the noise floor and find the intersection
%            time T_i where the decaying signal meets the noise floor.
%          • Truncate the backward integral at T_i.
%          • Add the Chu correction C = noise_power × T_late / (6 ln10) × fs
%            to compensate for the missing exponential tail.
%          • T20 is then fitted over −5 to −25 dB on the corrected EDC and
%            extrapolated to 60 dB.  Returns NaN if dynamic range < 25 dB.
%        This matches the default 'cutWithCorrection' method of the ITA
%        Toolbox (ita_roomacoustics_EDC).  Without this compensation,
%        background noise inflates the Schroeder integral, biasing T20
%        high when SNR is limited.
%     4. Band-wise clarity C80 = 10·log10(E_early / E_late), where the
%        early/late boundary is TEARLY seconds (typically 0.08 s).
%     5. Wideband centre time Ts computed from the aligned, full-band RIR.
%     6. Timbre ratios from band decay times:
%          BR = (T125 + T250) / (T500 + T1000)   (bass ratio)
%          TR = (T2000 + T4000) / (T500 + T1000)  (treble ratio)
%
%   INPUTS
%     h                   : mono impulse response (vector), any length
%     fs                  : sampling rate [Hz]
%     bandsPerOctave      : 1 for octave bands, 3 for third-octave bands
%     freqLims            : [fMin fMax] analysis band limits [Hz]
%     tEarly              : early/late boundary for C80 [s]  (typically 0.08)
%     alignThreshold_dB   : onset detection level below peak [dB]  (e.g. -20)
%
%   OUTPUTS
%     results : struct with fields
%       .bandCentersHz  — analysis band centre frequencies [Hz]
%       .T20            — decay time per band [s] (T20 estimator, ISO 3382)
%       .C80            — clarity per band [dB]
%       .Ts             — wideband centre time [s]
%       .BR             — bass ratio (NaN if required bands not present)
%       .TR             — treble ratio (NaN if required bands not present)
%       .hAligned       — onset-aligned impulse response (column vector)
%
%   DEPENDENCIES
%     Signal Processing Toolbox  (designfilt, filtfilt)
%
%   NOTES
%     • Bandpass filters are 8th-order Butterworth IIR, applied zero-phase
%       via filtfilt.  Bands whose edges violate the Nyquist limit are skipped.
%     • BR and TR require the relevant standard octave-band centres to fall
%       within FREQLIMS.  If any band is missing or its decay time is NaN,
%       the corresponding ratio is returned as NaN.
%     • The alignment pre-offset (5 samples) is fixed; adjust the constant
%       ALIGN_PRE_SAMPLES below if needed.
%
%   See also: polymaxRIR, polymax

% =========================================================================
%  Constants
% =========================================================================

ALIGN_PRE_SAMPLES = 5;    % samples prepended before detected onset

% =========================================================================
%  Input validation and setup
% =========================================================================

h    = h(:);
fMin = freqLims(1);
fMax = freqLims(2);

% =========================================================================
%  1. Onset alignment
% =========================================================================

hAligned = alignOnset(h, alignThreshold_dB, ALIGN_PRE_SAMPLES);
N = numel(hAligned);
t = (0 : N-1).' / fs;

% =========================================================================
%  2. Crop to mid-frequency T60 length
% =========================================================================
%
%   Samples beyond one T60 carry negligible energy and only increase
%   processing time.  The mid-frequency T60 is estimated from a quick
%   T20 Schroeder integration over a representative mid-frequency
%   range (500–1000 Hz).  The signal is then hard-trimmed to that
%   length before the full band-wise analysis runs.

MID_LO = 500;   % [Hz]  lower edge of mid-frequency range
MID_HI = 1000;  % [Hz]  upper edge of mid-frequency range

T60mid_crop = NaN;
if MID_HI < fs/2
    try
        bpMid = designfilt('bandpassiir', ...
            'FilterOrder',         8, ...
            'HalfPowerFrequency1', MID_LO, ...
            'HalfPowerFrequency2', MID_HI, ...
            'SampleRate',          fs);
        hMid      = filtfilt(bpMid, hAligned);
        T60mid_crop = simpleT20(hMid.^2, (0:numel(hMid)-1).'/fs, fs);
    catch
        % If the quick estimate fails, skip cropping and proceed in full.
    end
end

if ~isnan(T60mid_crop) && T60mid_crop > 0
    NcropMax = round(T60mid_crop * fs);
    if NcropMax < numel(hAligned)
        fprintf('acoustic_indices: cropping IR %d -> %d samples  (T60_mid = %.3f s)\n', ...
            numel(hAligned), NcropMax, T60mid_crop);
        hAligned = hAligned(1:NcropMax);
    end
else
    fprintf('acoustic_indices: mid-T60 estimate unavailable — no cropping applied\n');
end

% Recompute N and t after potential crop
N = numel(hAligned);
t = (0 : N-1).' / fs;

% =========================================================================
%  3. Analysis band definitions
% =========================================================================

%   Centre frequencies follow the exact ISO 266 / IEC 61260 progression:
%     fc = 1000 * 2^(n/B),  with edges  fc / 2^(1/(2B))  and  fc * 2^(1/(2B)).
%   Bands are selected by centre frequency falling within [fMin, fMax].
%   Edges are NOT clipped to fMin/fMax so that each band has standard ISO width.
fRef     = 1000;
halfStep = 2^(1 / (2*bandsPerOctave));
nMin     = floor(bandsPerOctave * log2(fMin / fRef)) - 1;
nMax     = ceil( bandsPerOctave * log2(fMax / fRef)) + 1;
n_all    = nMin : nMax;
fc_all   = fRef * 2.^(n_all / bandsPerOctave);
fc       = fc_all(fc_all >= fMin & fc_all <= fMax);   % centre-based selection
nBands   = numel(fc);

T20 = nan(nBands, 1);
C80 = nan(nBands, 1);

% =========================================================================
%  4. Wideband centre time Ts
% =========================================================================

eWide = hAligned.^2;
if sum(eWide) > 0
    Ts = sum(t .* eWide) / sum(eWide);
else
    Ts = NaN;
end

% =========================================================================
%  5. Band-wise decay time and clarity
% =========================================================================

for i = 1 : nBands
    fci = fc(i);
    % ISO band edges: exact geometric means between consecutive centres
    f1 = fci / halfStep;
    f2 = fci * halfStep;

    % Skip bands that violate Nyquist or zero-frequency
    if f1 <= 0 || f2 >= fs/2
        continue;
    end

    % 8th-order zero-phase Butterworth bandpass
    bpFilt = designfilt('bandpassiir', ...
        'FilterOrder',         8, ...
        'HalfPowerFrequency1', f1, ...
        'HalfPowerFrequency2', f2, ...
        'SampleRate',          fs);

    hBand = filtfilt(bpFilt, hAligned);
    eBand = hBand.^2;

    % T20 from Schroeder integration
    T20(i) = estimateDecayTime(eBand, fs);

    % Clarity: early-to-late energy ratio with boundary at tEarly
    C80(i) = computeClarity(eBand, fs, tEarly);
end

% =========================================================================
%  6. Timbre ratios  (BR, TR)
% =========================================================================
%
%   BR = (T125 + T250)  / (T500 + T1000)
%   TR = (T2000 + T4000) / (T500 + T1000)

T125  = getBandValue(fc, T20, 125,  bandsPerOctave);
T250  = getBandValue(fc, T20, 250,  bandsPerOctave);
T500  = getBandValue(fc, T20, 500,  bandsPerOctave);
T1000 = getBandValue(fc, T20, 1000, bandsPerOctave);
T2000 = getBandValue(fc, T20, 2000, bandsPerOctave);
T4000 = getBandValue(fc, T20, 4000, bandsPerOctave);

denom = T500 + T1000;

if any(isnan([T125, T250, T500, T1000])) || denom <= 0
    BR = NaN;
else
    BR = (T125 + T250) / denom;
end

if any(isnan([T2000, T4000, T500, T1000])) || denom <= 0
    TR = NaN;
else
    TR = (T2000 + T4000) / denom;
end

% =========================================================================
%  Pack results
% =========================================================================

results.bandCentersHz = fc(:);
results.T20           = T20;
results.C80           = C80;
results.Ts            = Ts;
results.BR            = BR;
results.TR            = TR;
results.hAligned      = hAligned;

% =========================================================================
%  Diagnostic plots
% =========================================================================

hOriginal = h;
t_orig    = (0 : numel(hOriginal)-1).' / fs;
t_aligned = (0 : numel(hAligned)-1).'  / fs;

figure('Name', 'acoustic_indices — aligned RIR', 'Color', 'w');
plot(t_orig,    hOriginal,  'DisplayName', 'Original h');  hold on;
plot(t_aligned, hAligned,   'DisplayName', 'Aligned h');
grid on;
xlabel('Time (s)');
ylabel('Amplitude');
title('Original vs onset-aligned impulse response');
legend('show', 'Location', 'best');

figure('Name', 'acoustic_indices — decay and clarity', 'Color', 'w');

subplot(2, 1, 1);
semilogx(fc, T20, '-o', 'LineWidth', 1.2, 'MarkerSize', 5);
grid on;
xlabel('Frequency (Hz)');
ylabel('Decay time (s)');
title('Reverberation time T_{20} per band  (ISO 3382 estimate of T_{60})');

subplot(2, 1, 2);
semilogx(fc, C80, '-o', 'LineWidth', 1.2, 'MarkerSize', 5);
grid on;
xlabel('Frequency (Hz)');
ylabel('C80 (dB)');
title(sprintf('C80 per band  (early boundary = %.0f ms)', 1000*tEarly));

fprintf('--- acoustic_indices results ---\n');
disp(results);

end   % acoustic_indices

% =========================================================================
%  Local helper: alignOnset
% =========================================================================

function hOut = alignOnset(h, threshold_dB, preSamples)
%ALIGNONSET  Trim leading silence before the direct-sound arrival.
%
%   Detects the first sample whose magnitude is at least THRESHOLD_DB below
%   the global peak, then shifts back by PRESAMPLES to preserve the onset
%   transient.  If detection fails, h is returned unchanged.

h   = h(:);
mag = abs(h);
mx  = max(mag);

if mx == 0
    hOut = h;
    return;
end

thr_lin = mx * 10^(threshold_dB / 20);   % linear threshold
idx     = find(mag >= thr_lin, 1, 'first');

if isempty(idx)
    hOut = h;
    return;
end

idx0 = max(1, idx - preSamples);
hOut = h(idx0 : end);
end

% =========================================================================
%  Local helper: estimateDecayTime
% =========================================================================

function T = estimateDecayTime(e, fs)
%ESTIMATEDECAYTIME  T20 from noise-compensated Schroeder back-integration.
%
%   Implements the Lundeby (1995) iterative noise-compensation algorithm
%   followed by the Chu (1978) correction term, matching the default
%   'cutWithCorrection' method of the ITA Toolbox.
%
%   Falls back to a plain Schroeder integral (simpleT20) when the signal
%   is too short or the Lundeby iteration cannot converge.

e = e(:);
N = numel(e);
if N < 10,  T = NaN;  return;  end
t = (0 : N-1).' / fs;

% ---- constants (matching ITA defaults) ----------------------------------
nPartsPer10dB    = 5;    % smoothing intervals per 10 dB of decay
dbAboveNoise     = 10;   % regression stops this many dB above noise floor
dynForReg        = 20;   % regression window width [dB]
MAX_ITER         = 30;

% ---- 1. Initial block-average (~30 ms blocks) ---------------------------
nBlk  = max(1, round(0.030 * fs));
nFull = floor(N / nBlk);
if nFull < 4
    T = simpleT20(e, t, fs);
    return;
end
[eBlk, tBlk] = makeBlocks(e, N, nBlk, fs);

% ---- 2. Initial noise estimate (last 10 % of blocks) --------------------
noiseEst = estimateNoise(eBlk);

% ---- 3. Preliminary regression ------------------------------------------
[~, iPeak] = max(eBlk);
iStop = find(10*log10(eBlk(iPeak+1:end) + realmin) > ...
             10*log10(noiseEst) + dbAboveNoise, 1, 'last');
if isempty(iStop) || iStop < 2
    T = simpleT20(e, t, fs);
    return;
end
iStop = iStop + iPeak;
iStop = min(iStop, numel(eBlk));   % clamp

[c, ok] = regressBlocks(eBlk, tBlk, iPeak, iStop);
if ~ok
    T = simpleT20(e, t, fs);
    return;
end
crossPt = (10*log10(noiseEst) - c(1)) / c(2);

% ---- 4. Refine block size before iteration ------------------------------
%   nPartsPer10dB blocks per 10 dB of decay over [iPeak, iStop]
dynRange_prelim = abs(diff(10*log10(eBlk([iPeak, iStop]) + realmin)));
if dynRange_prelim < 1,  dynRange_prelim = 1;  end
nBlkInDecay = max(1, dynRange_prelim / 10 * nPartsPer10dB);
timeSpan    = tBlk(iStop) - tBlk(iPeak);
nBlk        = max(1, min(round(timeSpan / nBlkInDecay * fs), round(0.1*fs)));
nFull       = floor(N / nBlk);
if nFull < 4
    T = simpleT20(e, t, fs);
    return;
end
[eBlk, tBlk] = makeBlocks(e, N, nBlk, fs);
[~, iPeak]   = max(eBlk);

% ---- 5. Lundeby iteration -----------------------------------------------
for iter = 1 : MAX_ITER

    % Estimate noise: from last 10% or 10 dB below crossing point
    idx10dBBefore = max(1, round((crossPt - dynForReg / max(abs(c(2)), 0.1)) * fs / nBlk));
    idxNoise90    = max(1, round(0.9 * numel(eBlk)));
    idxNoiseStart = min(idxNoise90, idx10dBBefore);
    noiseEst      = estimateNoise(eBlk(idxNoiseStart:end));

    % Find regression window
    iStart2 = find(10*log10(eBlk(iPeak:end) + realmin) < ...
                   10*log10(noiseEst) + dbAboveNoise + dynForReg, 1, 'first');
    if isempty(iStart2),  break;  end
    iStart2 = min(iStart2 + iPeak - 1, numel(eBlk));

    iStop2 = find(10*log10(eBlk(iStart2+1:end) + realmin) < ...
                  10*log10(noiseEst) + dbAboveNoise, 1, 'first');
    if isempty(iStop2),  break;  end
    iStop2 = min(iStop2 + iStart2, numel(eBlk));   % clamp

    if iStop2 <= iStart2,  break;  end

    [c, ok] = regressBlocks(eBlk, tBlk, iStart2, iStop2);
    if ~ok,  break;  end

    oldCrossPt = crossPt;
    crossPt    = (10*log10(noiseEst) - c(1)) / c(2);

    if abs(crossPt - oldCrossPt) < 0.01,  break;  end
end

% ---- 6. Truncate and apply Chu correction -------------------------------
crossPt = max(0, min(crossPt, t(end)));
iCut    = min(N, max(1, round(crossPt * fs)));

T_late  = -60 / min(c(2), -0.01);          % T60 from regression slope
C_chu   = noiseEst * T_late / (6*log(10)) * fs;

% ---- 7. Corrected EDC and T20 fit ---------------------------------------
E  = flipud(cumsum(e(iCut:-1:1))) + C_chu;
mx = max(E);
if mx <= 0,  T = NaN;  return;  end
EdB = 10*log10(E / mx + eps);
tCut = t(1:iCut);

tailIdx = max(1, numel(EdB) - floor(0.3*numel(EdB)) + 1) : numel(EdB);
if -median(EdB(tailIdx)) < 25,  T = NaN;  return;  end

T = fitSchroederSlope(EdB, tCut, -5, -25, 60);
end

% ---- helpers ------------------------------------------------------------

function [eB, tB] = makeBlocks(e, N, nBlk, fs)
nF  = floor(N / nBlk);
eB  = sum(reshape(e(1:nF*nBlk), nBlk, nF), 1).' / nBlk;
tB  = ((0:nF-1).' + 0.5) * nBlk / fs;
end

function n = estimateNoise(eBlk)
n = max(mean(eBlk(max(1, round(0.9*numel(eBlk))):end)), realmin);
end

function [c, ok] = regressBlocks(eBlk, tBlk, i1, i2)
i1 = max(1, i1);  i2 = min(numel(eBlk), i2);
ok = false;  c = [0; -1];
if i2 <= i1,  return;  end
X = [ones(i2-i1+1, 1), tBlk(i1:i2)];
c = X \ (10*log10(eBlk(i1:i2) + realmin));
ok = (c(2) < 0);
end

% -------------------------------------------------------------------------

function T = simpleT20(e, t, fs)
%SIMPLET20  Fallback: plain Schroeder integral without noise compensation.
%   Used when the signal is too short or the Lundeby iteration degenerates.
E = flipud(cumsum(flipud(e)));
mx = max(E);
if mx <= 0, T = NaN; return; end
EdB = 10*log10(E / mx + eps);
tailIdx = max(1, numel(EdB) - floor(0.3*numel(EdB)) + 1) : numel(EdB);
if -median(EdB(tailIdx)) < 25, T = NaN; return; end
T = fitSchroederSlope(EdB, t, -5, -25, 60);
end

% -------------------------------------------------------------------------

function T = fitSchroederSlope(EdB, t, upper_dB, lower_dB, delta_dB)
%FITSCHROEDERSLOPE  Linear fit to a dB segment of the Schroeder curve.
%
%   Fits a straight line to the segment between UPPER_DB and LOWER_DB and
%   extrapolates the slope to a total drop of DELTA_DB dB.

idx = find(EdB <= upper_dB & EdB >= lower_dB);
if numel(idx) < 2
    T = NaN;
    return;
end
p     = polyfit(t(idx), EdB(idx), 1);
slope = p(1);
if slope < 0
    T = -delta_dB / slope;
    if ~isfinite(T) || T <= 0, T = NaN; end
else
    T = NaN;
end
end

% =========================================================================
%  Local helper: computeClarity
% =========================================================================

function C = computeClarity(e, fs, tBoundary)
%COMPUTECLARITY  Early-to-late energy ratio in dB.
%
%   C = 10 * log10( E_early / E_late )
%   where the early window is [0, tBoundary] and the late window is
%   (tBoundary, end].  Returns +50 dB when late energy is negligible.

e = e(:);
t = (0 : numel(e)-1).' / fs;

E_early = sum(e(t <= tBoundary));
E_late  = sum(e(t >  tBoundary));

if E_late <= 0
    C = 50;   % effectively infinite clarity
    return;
end
C = 10*log10(E_early / E_late);
end

% =========================================================================
%  Local helper: getBandValue
% =========================================================================

function Tt = getBandValue(fc, Tvec, fTarget, bandsPerOctave)
%GETBANDVALUE  Return the decay-time entry matching a target band centre.
%
%   Searches FC for the closest entry to FTARGET in log-frequency, within
%   a tolerance of half a band step (with a small cushion).  Returns NaN if
%   no band matches (i.e., FTARGET is outside FREQLIMS).

fc   = fc(:);
Tvec = Tvec(:);

tol_log2 = (1 / (2*bandsPerOctave)) * 0.51;   % half-step + cushion

[dmin, idx] = min(abs(log2(fc / fTarget)));

if isempty(idx) || dmin > tol_log2
    Tt = NaN;
    return;
end

Tt = Tvec(idx);
if ~isfinite(Tt)
    Tt = NaN;
end
end