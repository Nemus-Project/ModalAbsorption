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
%     2. Band-wise decay time from Schroeder back-integration:
%          • T30 (fit over -5 to -35 dB) when dynamic range >= 35 dB.
%          • T20 (fit over -5 to -25 dB) as fallback when range >= 25 dB.
%          • NaN otherwise.
%        All values are stored in RESULTS.T30 regardless of which estimator
%        was used; the adaptive choice is transparent to the caller.
%     3. Band-wise clarity C80 = 10·log10(E_early / E_late), where the
%        early/late boundary is TEARLY seconds (typically 0.08 s).
%     4. Wideband centre time Ts computed from the aligned, full-band RIR.
%     5. Timbre ratios from band decay times:
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
%       .T30            — decay time per band [s] (T30 or T20 adaptively)
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
%  2. Analysis band definitions
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

T30 = nan(nBands, 1);
C80 = nan(nBands, 1);

% =========================================================================
%  3. Wideband centre time Ts
% =========================================================================

eWide = hAligned.^2;
if sum(eWide) > 0
    Ts = sum(t .* eWide) / sum(eWide);
else
    Ts = NaN;
end

% =========================================================================
%  4. Band-wise decay time and clarity
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

    % Adaptive T30/T20 from Schroeder integration
    T30(i) = estimateDecayTime(eBand, fs);

    % Clarity: early-to-late energy ratio with boundary at tEarly
    C80(i) = computeClarity(eBand, fs, tEarly);
end

% =========================================================================
%  5. Timbre ratios  (BR, TR)
% =========================================================================
%
%   BR = (T125 + T250)  / (T500 + T1000)
%   TR = (T2000 + T4000) / (T500 + T1000)

T125  = getBandValue(fc, T30, 125,  bandsPerOctave);
T250  = getBandValue(fc, T30, 250,  bandsPerOctave);
T500  = getBandValue(fc, T30, 500,  bandsPerOctave);
T1000 = getBandValue(fc, T30, 1000, bandsPerOctave);
T2000 = getBandValue(fc, T30, 2000, bandsPerOctave);
T4000 = getBandValue(fc, T30, 4000, bandsPerOctave);

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
results.T30           = T30;
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
semilogx(fc, T30, '-o', 'LineWidth', 1.2, 'MarkerSize', 5);
grid on;
xlabel('Frequency (Hz)');
ylabel('Decay time (s)');
title('Per-band decay time (T30 preferred, T20 fallback)');

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
%ESTIMATEDECAYTIME  Adaptive T30/T20 from Schroeder back-integration.
%
%   Integrates the squared band-filtered impulse response backward in time
%   (Schroeder integral), normalises to 0 dB, then estimates the noise floor
%   from the last 30 % of the record to determine usable dynamic range.
%
%   Decision logic:
%     dynRange >= 35 dB  ->  T30  (fit -5 to -35 dB, extrapolated x2)
%     dynRange >= 25 dB  ->  T20  (fit -5 to -25 dB, extrapolated x3)
%     otherwise          ->  NaN

e = e(:);
N = numel(e);
if N < 10
    T = NaN;
    return;
end

t = (0 : N-1).' / fs;

% Schroeder integral (backward cumulative sum), normalised
E   = flipud(cumsum(flipud(e)));
mx  = max(E);
if mx <= 0
    T = NaN;
    return;
end
EdB = 10*log10(E / mx + eps);   % peak ≈ 0 dB

% Estimate noise floor from the last 30 % of the record
tailIdx      = max(1, N - floor(0.3*N) + 1) : N;
noiseFloor   = median(EdB(tailIdx));
dynRange_dB  = -noiseFloor;      % positive value

T = NaN;
if dynRange_dB >= 35
    T = fitSchroederSlope(EdB, t, -5, -35, 30);
    if ~isnan(T), return; end
end
if dynRange_dB >= 25
    T = fitSchroederSlope(EdB, t, -5, -25, 20);
end
end

% -------------------------------------------------------------------------

function T = fitSchroederSlope(EdB, t, upper_dB, lower_dB, delta_dB)
%FITSCHROEDERSLOPE  Linear fit to a dB segment of the Schroeder curve.
%
%   Fits a straight line to the segment [upper_dB, lower_dB] and extrapolates
%   the slope to a DELTA_DB drop to recover the decay time.

idx = find(EdB <= upper_dB & EdB >= lower_dB);
if numel(idx) < 2
    T = NaN;
    return;
end
p     = polyfit(t(idx), EdB(idx), 1);   % p(1) = slope [dB/s], negative
slope = p(1);
if slope < 0
    T = -delta_dB / slope;
    if ~isfinite(T) || T <= 0
        T = NaN;
    end
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