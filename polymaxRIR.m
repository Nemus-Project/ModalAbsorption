%==========================================================================
% polymaxRIR.m
%==========================================================================
% PolyMax-based iterative modal fitting for room impulse responses.
%
% TWO-REGIME STRATEGY
%   f < fTransition  (below Schroeder frequency)  — PolyMax modal fitting.
%     Polynomial order driven by physics (up to 6 × expected mode count),
%     capped by Nbins/4 and maxPolyOrderLF.  Poles returned directly from
%     the stability diagram (NoCluster = 1).
%
%   f >= fTransition (above Schroeder frequency)  — spectral peak-picking.
%     Modal overlap above Schroeder makes individual pole resolution
%     unreliable.  Peaks are detected in |H(f)| and assigned a uniform
%     damping derived from the Schroeder T30.
%
% WORKFLOW
%   1.  Load and normalise the .wav file.
%   2.  Compute acoustic indices (onset alignment, T30, C80, BR, TR).
%   3.  Compute the half-spectrum via FFT.
%   4.  Auto-estimate the Schroeder transition frequency if not set.
%   5.  Build two independent band grids:
%         • Processing bands — narrow linear bands of fixed Hz width.
%           PolyMax / peak-picking runs on these; small bandwidths keep
%           polynomial orders tractable regardless of bandMode.
%         • Display bands — exact ISO 266 octave or third-octave bands.
%           All statistics, console output, and figures use these.
%   6.  Run PolyMax (LF) or peak-picking (HF) per processing band.
%   7.  Least-squares amplitude solve; gain-outlier rejection.
%   8.  Optional iterative residual refinement.
%   9.  Optional FIR correction filter.
%  10.  Remap poles from processing bands to display bands; compute stats.
%  11.  Print summary and produce plots (in display-band resolution).
%  12.  Save results and normalised synthetic RIR.
%
% DEPENDENCIES
%   Signal Processing Toolbox        (filtfilt, designfilt, findpeaks)
%   Statistics and Machine Learning  (linkage, cluster)
%   Audio Toolbox                    (audioread / audiowrite)
%   polymax.m, acoustic_indices.m    (companion files in the same folder)
%
% Author:  Michele Ducceschi, March 2026
%==========================================================================

clear; close all; clc;

%% ========================== USER PARAMETERS ==============================

% --- Input folder ---------------------------------------------------------
%   INPUTFOLDER is the directory containing all .wav files to process.
%   Every .wav file found directly inside this folder will be analysed in
%   sequence.  The folder's base name is used as a subfolder within
%   outRootDir so that results from different measurement sets (empty room,
%   with specimen, etc.) are kept strictly separate.
%
%   Example layout:
%     ./measured_IRs/empty/    → results in outRootDir/empty/
%     ./measured_IRs/St Gobain/ → results in outRootDir/St Gobain/
inputFolder = './measured_IRs/empty';

% --- Input normalisation --------------------------------------------------
normaliseInput = true;

% --- Analysis frequency range ---------------------------------------------
fLow  = 20;
fHigh = 12000;

% --- Band subdivision — DISPLAY -------------------------------------------
%   Controls the ISO band resolution used for all statistics and plots.
%   'thirdOctave' | 'fullOctave' | 'custom'
bandMode        = 'thirdOctave';
customBandEdges = [50 100 200 400 800 1600 3200 4000];  % used only for 'custom'

% --- Band subdivision — PROCESSING ----------------------------------------
%   PolyMax runs on narrow linear bands of this fixed width [Hz].
%   Smaller values → lower polynomial orders → faster, more stable fits.
%   Must be wider than a few FFT bins; a value between 20–100 Hz works well.
procBandWidth_Hz = 25;   % [Hz]  linear processing band width

% --- Room geometry --------------------------------------------------------
V_room  = 277;   % [m^3]
c_sound = 343;   % [m/s]

% --- Transition frequency -------------------------------------------------
%   Below fTransition: PolyMax modal fitting.
%   Above fTransition: spectral peak-picking with Schroeder T30 damping.
%   Set to [] for automatic computation: 2000 * sqrt(T60_mid / V_room),
%   where T60_mid is the mean T60 in the 500–1000 Hz range from acoustic_indices.
fTransition = [];   % [Hz]

% --- PolyMax order cap (LF bands only) ------------------------------------
%   Order driven by mode count; Nbins/4 cap kicks in before maxPolyOrderLF
%   in most practical cases.
maxPolyOrderLF = 2000;   % hard ceiling

% --- PolyMax iteration count ----------------------------------------------
%   Step size = max(1, round((NordMax - NordMin) / maxPolyIters)).
maxPolyIters = 100;

% --- HF peak-picking (f >= fTransition) ----------------------------------
minPeakProm_dB = 0.5;   % [dB]  minimum peak prominence
minPeakDist_Hz = 0.5;   % [Hz]  minimum inter-peak spacing

% --- PolyMax stability tolerances -----------------------------------------
tolf_init   = 0.05;   % relative frequency tolerance, initial pass
told_init   = 0.10;   % relative damping  tolerance, initial pass
tolf_refine = 0.05;   % relative frequency tolerance, refinement passes
told_refine = 0.10;   % relative damping  tolerance, refinement passes

% --- T60 rejection ceiling ------------------------------------------------
T60cutMult  = 1.5;    % T60cut = T60cutMult × T60est  (per band)
T60cutFloor = 5.0;    % [s]  minimum allowable T60cut

% --- Live stabilisation diagram -------------------------------------------
%   true = draw one figure per band while PolyMax runs (debug / inspection)
LivePlot = false;

% --- Least-squares amplitude solve ----------------------------------------
overlapFrac = 0.10;   % fractional band overlap for LS window

% --- Gain outlier rejection -----------------------------------------------
enableGainOutlierRejection = true;
gainOutlierMethod          = 'mad';    % 'mad' | 'std'
gainOutlierZthresh         = 3.0;
maxGainOutlierLSIters      = 5;
minModesToKeepAfterReject  = 3;

% --- Iterative residual refinement ----------------------------------------
useIterativeRefinement = false;
nRefinementIters       = 8;
mergeTol_Hz            = 0.25;    % [Hz]  poles closer than this are merged
minNRMSEImprovement    = 1e-5;
maxNoImproveIters      = 2;

% --- FIR correction filter ------------------------------------------------
useFIRCorrection = false;
FIR_Ntaps        = 256;
FIR_fitSeconds   = 0.15;    % [s]  time window used for FIR LS fit
FIR_lambda       = 1e-6;    % Tikhonov regularisation for FIR solve

% --- Acoustic indices -----------------------------------------------------
MinSecondsForFFT  = 5;       % minimum FFT length [s] (zero-padded if needed)
bandsPerOctave    = 3;
freqLimsAI        = [20, fHigh];
tEarly            = 0.080;    % [s]  C80 early/late boundary
alignThreshold_dB = -5;       % [dB]  onset detection level below peak

% --- Quality threshold (keep vs discard decision) -------------------------
nRMSE_time_keep_threshold = 5e-3;

% --- Output ---------------------------------------------------------------
verbose    = true;
outRootDir = './polymaxRIR_results/';

%% ========================== LOAD AND PREPARE ==============================
%
%   Discover all .wav files in INPUTFOLDER, create the output subfolder
%   named after the input folder, then iterate over every file.  Each file
%   is processed independently; errors are caught per-file so a single
%   bad recording does not abort the whole batch.

% --- Discover .wav files --------------------------------------------------
wavList = dir(fullfile(inputFolder, '*.wav'));
if isempty(wavList)
    error('polymaxRIR:noWavFiles', ...
        'No .wav files found in ''%s''.', inputFolder);
end
NFiles = numel(wavList);
fprintf('\n=== polymaxRIR  BATCH MODE ===\n');
fprintf('Input folder : %s\n', inputFolder);
fprintf('Files found  : %d\n\n', NFiles);

% --- Output subfolder named after the input folder ------------------------
[~, folderName] = fileparts(strtrim(inputFolder));
if isempty(folderName)   % handle trailing path separator
    [~, folderName] = fileparts(fileparts(inputFolder));
end
outSubDir = fullfile(outRootDir, folderName);
if ~exist(outSubDir, 'dir'), mkdir(outSubDir); end

% --- Save the user-supplied fTransition so it can be reset each iteration -
fTransition_user = fTransition;

% --- Batch summary initialisation -----------------------------------------
%   NOTE: do NOT use struct('field', cell(N,1), ...) here — that syntax
%   creates an N-element struct array, not a scalar struct with array fields.
%   Field-by-field assignment produces the correct scalar struct.
batchSummary.filename       = cell(NFiles, 1);
batchSummary.nRMSE_IIR      = nan(NFiles, 1);
batchSummary.nRMSE_final    = nan(NFiles, 1);
batchSummary.totalPoles     = nan(NFiles, 1);
batchSummary.fTransition    = nan(NFiles, 1);
batchSummary.T60global_mean = nan(NFiles, 1);
batchSummary.isKept         = false(NFiles, 1);
batchSummary.errorMsg       = cell(NFiles, 1);

%% ========================= PER-FILE PROCESSING LOOP ======================

for iFile = 1 : NFiles

filename    = wavList(iFile).name(1 : end-4);   % strip .wav extension
fTransition = fTransition_user;                  % reset for each file

fprintf('\n--- File %d/%d : %s ---\n', iFile, NFiles, filename);

try   % catch errors per-file so the batch continues on failure

filenamewav = fullfile(inputFolder, wavList(iFile).name);
[inputRaw, fs] = audioread(filenamewav);

% Use first channel only; remove DC offset
if size(inputRaw, 2) > 1, inputRaw = inputRaw(:, 1); end
inputRaw = inputRaw - mean(inputRaw);

if normaliseInput
    pk = max(abs(inputRaw));
    if pk > 0, inputRaw = inputRaw / pk; end
    fprintf('Input normalised  (peak = %.4g)\n', pk);
end

fprintf('fs   : %d Hz  |  N = %d samples (%.2f s)\n', ...
    fs, numel(inputRaw), numel(inputRaw)/fs);

%% ===================== ACOUSTIC INDICES AND ALIGNED IR ====================

refAI = acoustic_indices(inputRaw, fs, bandsPerOctave, ...
    freqLimsAI, tEarly, alignThreshold_dB);

fprintf('BR (ref) = %.4g  |  TR (ref) = %.4g\n', refAI.BR, refAI.TR);

% Work on the onset-aligned IR from here on
input = refAI.hAligned(:);

% Build per-band damping look-up from T30 estimates (used later for T60cut
% and HF peak-picking).  T20 is an ISO 3382 estimate of T60 (decay time
% fitted over -5 to -35 dB and extrapolated to 60 dB), so T60 = T30 here.
fcBandsAI  = refAI.bandCentersHz(:);
T20bandsAI = refAI.T20(:);
idxV       = ~isnan(T20bandsAI) & (T20bandsAI > 0);
T60v       = T20bandsAI(idxV);   % T20 is already an estimate of T60
cBandsAI   = nan(size(T20bandsAI));
cBandsAI(idxV) = 3*log(10) ./ T60v;

% Fallback damping coefficient for bands where T30 is unavailable
cFallback = mean(cBandsAI(idxV), 'omitnan');
if isnan(cFallback) || cFallback <= 0, cFallback = 1; end

fprintf('Aligned IR : %d samples (%.3f s)  |  cFallback = %.4g\n', ...
    numel(input), numel(input)/fs, cFallback);

%% ============================= FFT ========================================

Lmin = round(MinSecondsForFFT * fs);
Lfft = max(Lmin, numel(input));
[spcTarget, fv, ~] = takefft(input, fs, Lfft);
omvFull = 2*pi*fv;
df = fv(2) - fv(1);

fprintf('FFT : Lfft = %d  |  bins = %d  |  df = %.4f Hz\n', ...
    Lfft, numel(spcTarget), df);

%% ========================== TRANSITION FREQUENCY ==========================

if isempty(fTransition)
    % Use the mid-frequency T60 (500–1000 Hz) from acoustic_indices.
    % This is the standard room-acoustics convention for the Schroeder frequency.
    fcValid   = fcBandsAI(idxV);          % centres corresponding to T60v entries
    idxMid    = (fcValid >= 500) & (fcValid <= 1000);
    T60mid    = mean(T60v(idxMid), 'omitnan');
    if isnan(T60mid) || T60mid <= 0
        error('polymaxRIR:noT60mid', ...
            ['acoustic_indices returned no valid T30 values in the 500–1000 Hz ' ...
             'range — cannot compute the Schroeder frequency.  Either supply ' ...
             'fTransition manually or widen freqLimsAI to cover 500–1000 Hz.']);
    end
    fTransition = 2000 * sqrt(T60mid / V_room);
    fprintf('fTransition (Schroeder, from acoustic_indices 500–1000 Hz) : %.1f Hz  (T60_mid = %.2f s)\n', ...
        fTransition, T60mid);
else
    fprintf('fTransition (user supplied) : %.1f Hz\n', fTransition);
end

%% ========================== BAND DEFINITION ===============================

% --- Display bands (ISO octave / third-octave) ----------------------------
%   Used for all statistics, console output, and figures.
switch lower(bandMode)
    case 'thirdoctave'
        [dispEdgesVec, dispCentresISO] = buildISOBands(fLow, fHigh, 3);
    case 'fulloctave'
        [dispEdgesVec, dispCentresISO] = buildISOBands(fLow, fHigh, 1);
    case 'custom'
        dispEdgesVec   = customBandEdges(:).';
        dispEdgesVec   = dispEdgesVec(dispEdgesVec >= fLow & dispEdgesVec <= fHigh);
        if numel(dispEdgesVec) < 2
            error('polymaxRIR:customBands', ...
                'customBandEdges must have >= 2 values within [fLow, fHigh].');
        end
        dispCentresISO = sqrt(dispEdgesVec(1:end-1) .* dispEdgesVec(2:end));
    otherwise
        error('polymaxRIR:unknownBandMode', 'Unknown bandMode ''%s''.', bandMode);
end

% Clip display bands to the usable FFT range (drop whole bands, not edges)
dispEdgesVec = dispEdgesVec(:).';
keepDisp     = (dispEdgesVec(2:end) > fv(2)) & (dispEdgesVec(1:end-1) < fv(end));
if ~all(keepDisp)
    dIdx         = find(keepDisp);
    dispEdgesVec = [dispEdgesVec(dIdx), dispEdgesVec(dIdx(end) + 1)];
end
% Recompute centres from surviving edges (guards against any edge dropping)
dispCentresISO = sqrt(dispEdgesVec(1:end-1) .* dispEdgesVec(2:end));
NDispBands     = numel(dispEdgesVec) - 1;

fprintf('Display bands : mode = %s  |  N = %d  |  [%.2f – %.2f] Hz\n', ...
    bandMode, NDispBands, dispEdgesVec(1), dispEdgesVec(end));

% --- Processing bands (linear, fixed width) --------------------------------
%   PolyMax and peak-picking run on these narrow bands.  Their width is
%   procBandWidth_Hz, independent of bandMode, so polynomial orders remain
%   tractable even for wide display octave bands.
procEdgesVec = buildLinearBands(fLow, fHigh, procBandWidth_Hz, fv);
% Drop degenerate processing bands
procEdgesVec = procEdgesVec(:).';
validP       = find(diff(procEdgesVec) > 0.5*df);
if numel(validP) < numel(procEdgesVec) - 1
    procEdgesVec = [procEdgesVec(validP), procEdgesVec(validP(end) + 1)];
end
procCentres  = 0.5*(procEdgesVec(1:end-1) + procEdgesVec(2:end));
NProcBands   = numel(procEdgesVec) - 1;

fprintf('Processing bands : width = %.0f Hz  |  N = %d  |  [%.2f – %.2f] Hz\n', ...
    procBandWidth_Hz, NProcBands, procEdgesVec(1), procEdgesVec(end));

%% ===================== PRE-ALLOCATE STORAGE ==============================

MAX_POLES  = 5e5;
fAll       = zeros(MAX_POLES, 1);
cAll       = zeros(MAX_POLES, 1);
aAll       = zeros(MAX_POLES, 1);
bAll       = zeros(MAX_POLES, 1);
bandIdxAll = zeros(MAX_POLES, 1, 'int32');   % indexes into procEdgesVec
cumPoles   = 0;

% Per-processing-band struct (used during the fitting loop and refinement).
% Display-band statistics are computed separately after all poles are found.
poleStats_template = struct( ...
    'flow', NaN, 'fhigh', NaN, 'fcBand', NaN, 'N', 0, ...
    'densityRatio', NaN, 'modalDensity_Hz', NaN, ...
    'meanF', NaN, 'stdF', NaN, ...
    'meanSpacing_Hz', NaN, 'stdSpacing_Hz', NaN, 'cvSpacing', NaN, ...
    'meanT60', NaN, 'stdT60', NaN, 'cvT60', NaN, ...
    'meanC', NaN, 'stdC', NaN, 'regime', 'LF');
procPoleStats = repmat(poleStats_template, 1, NProcBands);

%% =================== INITIAL POLYMAX PASS (per processing band) ==========

fprintf('\n=== INITIAL POLYMAX PASS (%d processing bands) ===\n', NProcBands);

for ib = 1 : NProcBands
    flow   = procEdgesVec(ib);
    fhigh  = procEdgesVec(ib + 1);
    fcBand = procCentres(ib);
    bwBand = fhigh - flow;
    T60est = interpT60band(fcBand, fcBandsAI, T20bandsAI, cFallback);

    % Skip bands that are too narrow to support even a minimal fit
    minBW_Hz = max(3*df, 2.0);
    if bwBand < minBW_Hz
        if verbose
            fprintf('  Band %2d/%2d  [%6.1f – %6.1f Hz]  SKIPPED (BW = %.2f Hz)\n', ...
                ib, NProcBands, flow, fhigh, bwBand);
        end
        procPoleStats(ib) = emptyBandStats(flow, fhigh, fcBand, 'LF');
        continue;
    end

    % Extract the band's FFT slice
    [~, ibs] = min(abs(fv - flow));
    [~, ibe] = min(abs(fv - fhigh));
    ibe     = max(ibe, ibs + 1);
    fvBand  = fv(ibs:ibe);
    fftBand = spcTarget(ibs:ibe);
    Nbins   = numel(fvBand);

    regime = selectRegime(fcBand, fTransition);

    if strcmp(regime, 'LF')
        [pmNordMin, pmNordMax, pmStep, pmT60cut, pmNexpected] = ...
            autoPolymaxParams(fcBand, bwBand, Nbins, V_room, c_sound, T60est, ...
            maxPolyOrderLF, maxPolyIters, T60cutMult, T60cutFloor);
        if verbose
            fprintf(['  Band %2d/%2d  [%6.1f – %6.1f Hz]  LF  |  ' ...
                'fc = %.1f Hz  BW = %.1f Hz  bins = %d  |  ' ...
                'NordMax = %d  step = %d  Nexp = %d  T60cut = %.1f s\n'], ...
                ib, NProcBands, flow, fhigh, fcBand, bwBand, Nbins, ...
                pmNordMax, pmStep, pmNexpected, pmT60cut);
        end
        [fBand, cBand] = polymaxWrapper(fftBand, fvBand, ...
            pmNordMin, pmNordMax, pmStep, tolf_init, told_init, ...
            pmT60cut, pmNexpected, LivePlot);

    else   % HF peak-picking
        cFromT60 = 3*log(10) / max(T60est, 0.01);
        if verbose
            fprintf(['  Band %2d/%2d  [%6.1f – %6.1f Hz]  HF  |  ' ...
                'fc = %.1f Hz  BW = %.1f Hz  bins = %d  |  ' ...
                'peak-pick  T60est = %.2f s\n'], ...
                ib, NProcBands, flow, fhigh, fcBand, bwBand, Nbins, T60est);
        end
        [fBand, cBand] = peakPickWrapper(fftBand, fvBand, df, ...
            minPeakProm_dB, minPeakDist_Hz, cFromT60);
    end

    if isempty(fBand)
        if verbose, fprintf('    -> 0 poles\n'); end
        procPoleStats(ib) = emptyBandStats(flow, fhigh, fcBand, regime);
        continue;
    end

    % Capture the per-pole T60 spread BEFORE homogenisation.
    % After cBand(:) = cMean all values are identical, so computeBandStats
    % would return stdT60 = 0.  The spread from the stability diagram is
    % the physically meaningful quantity for the LF comparison plot.
    T60perPole = 3*log(10) ./ max(cBand, 1e-10);
    T60std_raw = std(T60perPole);
    if isnan(T60std_raw), T60std_raw = 0; end

    % Replace individual per-pole damping estimates with the band mean.
    % This prevents a handful of poorly-conditioned poles from producing
    % widely varying residuals when the amplitudes are solved below.
    cMean    = mean(cBand);
    cBand(:) = cMean;

    omvBand = 2*pi*fvBand(:);
    Np      = numel(fBand);
    [aBand, bBand] = leastSqWeights(omvBand, Np, fBand(:), cBand(:), fftBand(:));

    if enableGainOutlierRejection && Np > minModesToKeepAfterReject
        [aBand, bBand, fBand, cBand] = rejectGainOutliers( ...
            aBand, bBand, fBand, cBand, omvBand, fftBand, ...
            gainOutlierMethod, gainOutlierZthresh, ...
            maxGainOutlierLSIters, minModesToKeepAfterReject);
        Np = numel(fBand);
        if Np > 0
            cMean    = mean(cBand);
            cBand(:) = cMean;
        end
    end

    if Np == 0
        procPoleStats(ib) = emptyBandStats(flow, fhigh, fcBand, regime);
        continue;
    end

    fAll(cumPoles + 1 : cumPoles + Np)       = fBand(:);
    cAll(cumPoles + 1 : cumPoles + Np)       = cBand(:);
    aAll(cumPoles + 1 : cumPoles + Np)       = aBand(:);
    bAll(cumPoles + 1 : cumPoles + Np)       = bBand(:);
    bandIdxAll(cumPoles + 1 : cumPoles + Np) = ib;
    cumPoles = cumPoles + Np;

    procPoleStats(ib) = computeBandStats(fBand, cBand, flow, fhigh, fcBand, regime);
    procPoleStats(ib).stdT60 = T60std_raw;   % restore spread from original per-pole values

    if verbose
        fprintf('    -> %d poles  |  cMean = %.4g (T60 = %.3f s)  |  density = %.4f modes/Hz\n', ...
            Np, cMean, 3*log(10)/cMean, procPoleStats(ib).modalDensity_Hz);
    end
end

% Trim pre-allocated arrays to the number of poles found
fAll       = fAll(1:cumPoles);
cAll       = cAll(1:cumPoles);
aAll       = aAll(1:cumPoles);
bAll       = bAll(1:cumPoles);
bandIdxAll = bandIdxAll(1:cumPoles);

fprintf('\nInitial pass complete : %d total poles\n', cumPoles);

%% ===================== ITERATIVE RESIDUAL REFINEMENT =====================
%
%   outRef is simply the onset-aligned, normalised input signal.
%   There is no reason to round-trip through FFT/IFFT — that reintroduces
%   wrap-around artefacts for long-T60 signals.  The input is already the
%   ground truth we want to compare against.

outRef       = input(:);   % length = numel(input), no IFFT wrap-around
prevNRMSE    = Inf;
noImproveCnt = 0;

if useIterativeRefinement && nRefinementIters > 0
    fprintf('\n=== ITERATIVE RESIDUAL REFINEMENT ===\n');

    for iIter = 1 : nRefinementIters

        % Build synthetic spectrum band-by-band — used for the spectral residual only.
        spcSyn_cur = zeros(size(spcTarget));
        for ib = 1 : NProcBands
            idxIB = (bandIdxAll(1:cumPoles) == ib);
            if ~any(idxIB), continue; end
            [~, ibs_s] = min(abs(fv - procEdgesVec(ib)));
            [~, ibe_s] = min(abs(fv - procEdgesVec(ib + 1)));
            ibe_s = max(ibe_s, ibs_s + 1);
            omvSlice = 2*pi * fv(ibs_s:ibe_s);
            spcSyn_cur(ibs_s:ibe_s) = spcSyn_cur(ibs_s:ibe_s) + ...
                spectrumBuildSlice(aAll(idxIB), bAll(idxIB), ...
                                   cAll(idxIB), fAll(idxIB), omvSlice);
        end

        % nRMSE via direct time-domain synthesis — no IFFT wrap-around.
        outSyn_cur = synthesizeIR(aAll(1:cumPoles), bAll(1:cumPoles), ...
                                   cAll(1:cumPoles), fAll(1:cumPoles), fs, numel(outRef));
        curNRMSE   = normalisedRMSE_time(outRef, outSyn_cur);
        fprintf('\n  Iter %d/%d  |  nRMSE_time = %.6f\n', ...
            iIter, nRefinementIters, curNRMSE);

        % Early stopping: halt if improvement has stalled
        if (prevNRMSE - curNRMSE) < minNRMSEImprovement
            noImproveCnt = noImproveCnt + 1;
            fprintf('    Stalled (%d/%d)\n', noImproveCnt, maxNoImproveIters);
            if noImproveCnt >= maxNoImproveIters
                fprintf('    Stopping.\n');
                break;
            end
        else
            noImproveCnt = 0;
            prevNRMSE    = curNRMSE;
        end

        spcResidual = spcTarget - spcSyn_cur;
        nAdded      = 0;

        for ib = 1 : NProcBands
            flow   = procEdgesVec(ib);
            fhigh  = procEdgesVec(ib + 1);
            fcBand = procCentres(ib);
            bwBand = fhigh - flow;
            if bwBand < max(3*df, 2.0), continue; end

            T60est = interpT60band(fcBand, fcBandsAI, T20bandsAI, cFallback);

            [~, ibs] = min(abs(fv - flow));
            [~, ibe] = min(abs(fv - fhigh));
            ibe      = max(ibe, ibs + 1);
            fvBand   = fv(ibs:ibe);
            fftResid = spcResidual(ibs:ibe);
            Nbins    = numel(fvBand);

            regime = selectRegime(fcBand, fTransition);

            if strcmp(regime, 'LF')
                [pmNordMin, pmNordMax, pmStep, pmT60cut, pmNexpected] = ...
                    autoPolymaxParams(fcBand, bwBand, Nbins, V_room, c_sound, T60est, ...
                    maxPolyOrderLF, maxPolyIters, T60cutMult, T60cutFloor);
                [fNew, cNew] = polymaxWrapper(fftResid, fvBand, ...
                    pmNordMin, pmNordMax, pmStep, tolf_refine, told_refine, ...
                    pmT60cut, pmNexpected, false);
            else
                cFromT60 = 3*log(10) / max(T60est, 0.01);
                [fNew, cNew] = peakPickWrapper(fftResid, fvBand, df, ...
                    minPeakProm_dB, minPeakDist_Hz, cFromT60);
            end

            if isempty(fNew), continue; end

            % Reject new poles that are already represented in this band
            fExist     = fAll(bandIdxAll(1:cumPoles) == ib);
            [fNew, cNew] = mergeRejectClosePoles(fNew, cNew, fExist, mergeTol_Hz);
            if isempty(fNew), continue; end

            % Pull new damping toward existing band mean to stay consistent
            cExist   = cAll(bandIdxAll(1:cumPoles) == ib);
            cNew(:)  = mean([cExist(:); cNew(:)]);

            omvBand = 2*pi*fvBand(:);
            NpNew   = numel(fNew);
            [aNew, bNew] = leastSqWeights(omvBand, NpNew, fNew(:), cNew(:), fftResid(:));

            % Grow storage arrays if needed
            iEnd = cumPoles + NpNew;
            if iEnd > numel(fAll)
                fAll       = [fAll;       zeros(MAX_POLES, 1)];  %#ok<AGROW>
                cAll       = [cAll;       zeros(MAX_POLES, 1)];
                aAll       = [aAll;       zeros(MAX_POLES, 1)];
                bAll       = [bAll;       zeros(MAX_POLES, 1)];
                bandIdxAll = [bandIdxAll; zeros(MAX_POLES, 1, 'int32')];
            end

            fAll(cumPoles + 1 : iEnd)       = fNew(:);
            cAll(cumPoles + 1 : iEnd)       = cNew(:);
            aAll(cumPoles + 1 : iEnd)       = aNew(:);
            bAll(cumPoles + 1 : iEnd)       = bNew(:);
            bandIdxAll(cumPoles + 1 : iEnd) = ib;
            cumPoles = iEnd;
            nAdded   = nAdded + NpNew;
        end

        if nAdded == 0
            fprintf('    No new poles found. Stopping.\n');
            break;
        end
        fprintf('    Added %d new poles  (total : %d)\n', nAdded, cumPoles);

        % Pre-compute band membership once — reused by the damping re-average
        % and the stats update, replacing two O(NProcBands × cumPoles) scans.
        bandMembership = cell(NProcBands, 1);
        for ib = 1 : NProcBands
            bandMembership{ib} = find(bandIdxAll(1:cumPoles) == ib);
        end

        % Re-average damping per processing band after adding poles
        for ib = 1 : NProcBands
            idx = bandMembership{ib};
            if ~isempty(idx), cAll(idx) = mean(cAll(idx)); end
        end

        % Global banded LS re-solve for amplitudes (uses processing bands)
        [aAll(1:cumPoles), bAll(1:cumPoles)] = globalBandedLS( ...
            spcTarget, fv, fAll(1:cumPoles), cAll(1:cumPoles), ...
            bandIdxAll(1:cumPoles), procEdgesVec, overlapFrac);

        % Update per-processing-band statistics (regime info used in next iter)
        for ib = 1 : NProcBands
            idx = bandMembership{ib};
            if ~isempty(idx)
                procPoleStats(ib) = computeBandStats( ...
                    fAll(idx), cAll(idx), ...
                    procEdgesVec(ib), procEdgesVec(ib + 1), ...
                    procCentres(ib), procPoleStats(ib).regime);
            end
        end
    end
end

% Final trim and sort by frequency
fAll       = fAll(1:cumPoles);
cAll       = cAll(1:cumPoles);
aAll       = aAll(1:cumPoles);
bAll       = bAll(1:cumPoles);
bandIdxAll = bandIdxAll(1:cumPoles);
[fAll, iSort] = sort(fAll);
cAll = cAll(iSort); aAll = aAll(iSort);
bAll = bAll(iSort); bandIdxAll = bandIdxAll(iSort);

%% ========================= FINAL IIR SPECTRUM ============================
%
%   spcSyn_IIR is built band-by-band for the spectrum plot and residual.
%   The time-domain output uses direct modal synthesis — synthesizeIR —
%   so decays are physically correct regardless of FFT window length.

spcSyn_IIR = zeros(size(spcTarget));
for ib = 1 : NProcBands
    idxIB = (bandIdxAll == ib);
    if ~any(idxIB), continue; end
    [~, ibs_f] = min(abs(fv - procEdgesVec(ib)));
    [~, ibe_f] = min(abs(fv - procEdgesVec(ib + 1)));
    ibe_f = max(ibe_f, ibs_f + 1);
    omvSlice = 2*pi * fv(ibs_f:ibe_f);
    spcSyn_IIR(ibs_f:ibe_f) = spcSyn_IIR(ibs_f:ibe_f) + ...
        spectrumBuildSlice(aAll(idxIB), bAll(idxIB), ...
                           cAll(idxIB), fAll(idxIB), omvSlice);
end

Lout       = numel(outRef);
outSyn_IIR = synthesizeIR(aAll, bAll, cAll, fAll, fs, Lout);
nRMSE_IIR  = normalisedRMSE_time(outRef, outSyn_IIR);

fprintf('\nFinal IIR  |  poles = %d  |  nRMSE_time = %.6f\n', cumPoles, nRMSE_IIR);

%% ========================= FIR CORRECTION FILTER =========================

hFIR   = 1;      % identity filter (no FIR correction)
outSyn = outSyn_IIR;

if useFIRCorrection
    Nfit = min(Lout, round(FIR_fitSeconds * fs));
    if Nfit > FIR_Ntaps
        [hFIR, outSynFIR] = designResidualFIR_LS( ...
            outRef, outSyn_IIR, FIR_Ntaps, FIR_lambda, Nfit);
        outSyn = outSynFIR(1:Lout);
        fprintf('After FIR correction  |  nRMSE_time = %.6f\n', ...
            normalisedRMSE_time(outRef, outSyn));
    else
        warning('polymaxRIR:FIRskipped', ...
            'FIR correction skipped: Nfit (%d) <= FIR_Ntaps (%d).', ...
            Nfit, FIR_Ntaps);
    end
end

nRMSE_final = normalisedRMSE_time(outRef, outSyn);

%% ======= REMAP POLES TO DISPLAY BANDS AND COMPUTE DISPLAY STATISTICS =====
%
%   Now that all poles are finalised, assign each pole to its ISO display
%   band by frequency and compute per-display-band statistics.  These are
%   what is shown in all console output, figures, and the saved MAT file.

dispBandIdxAll = assignToDisplayBands(fAll, dispEdgesVec);

poleStats = computeDisplayBandStats( ...
    fAll, cAll, dispEdgesVec, dispCentresISO, dispBandIdxAll, fTransition);

fprintf('\n=== PER-BAND POLE STATISTICS (display bands) ===\n');
globalStats = computeGlobalStats(poleStats, fAll, cAll, dispEdgesVec, V_room, c_sound);
printAllStats(globalStats, poleStats, verbose);

%% ============================== PLOTS ====================================

tv = (0 : Lout - 1).' / fs;

% --- Spectrum and time-domain comparison ----------------------------------
hFig_spectrum = figure('Name', 'polymaxRIR — Spectrum', ...
    'Color', 'w', 'Position', [50 50 1100 700]);

subplot(3, 1, 1);
plot(fv, 20*log10(abs(spcTarget)  + eps), 'k', 'LineWidth', 0.7); hold on;
plot(fv, 20*log10(abs(spcSyn_IIR) + eps), 'r', 'LineWidth', 0.7);
for ib = 1 : NDispBands
    xline(dispEdgesVec(ib), '--', 'Color', [0.75 0.75 0.75], ...
        'HandleVisibility', 'off', 'LineWidth', 0.5);
end
xline(fTransition, '-', 'Color', [0 0.5 0], 'LineWidth', 1.2, 'Label', 'f_{Sch}');
xlim([fLow fHigh]);
xlabel('f (Hz)'); ylabel('|H| (dB)');
legend('Target', 'Synthetic', 'Location', 'best');
title(sprintf('%s — %d poles', filename, cumPoles));

subplot(3, 1, 2);
plot(fv, 20*log10(abs(spcTarget - spcSyn_IIR) + eps), ...
    'Color', [0 0.6 0], 'LineWidth', 0.7);
xlim([fLow fHigh]);
xlabel('f (Hz)'); ylabel('Residual (dB)');
title('Spectral residual');

subplot(3, 1, 3);
plot(tv*1e3, outRef,  'k', 'LineWidth', 0.6); hold on;
plot(tv*1e3, outSyn,  'r', 'LineWidth', 0.6);
xlabel('t (ms)'); ylabel('h(t)');
legend('Target', 'Synthetic', 'Location', 'best');
xlim([0, min(500, tv(end)*1e3)]);
title(sprintf('Time domain  (nRMSE = %.4g)', nRMSE_final));

% --- Per-band statistics --------------------------------------------------
hFig_bandStats = figure('Name', 'polymaxRIR — Per-band statistics', ...
    'Color', 'w', 'Position', [50 800 1200 600]);
plotBandStats(poleStats, globalStats, NDispBands, fTransition);

% --- LF T60 comparison: PolyMax poles vs acoustic_indices -----------------
hFig_LFcomp = figure('Name', 'polymaxRIR — LF T60 comparison', ...
    'Color', 'w', 'Position', [1300 800 800 500]);
plotLFComparison(poleStats, refAI, dispEdgesVec, dispCentresISO, fTransition);

%% ============================== SAVE =====================================

if ~exist(outSubDir, 'dir'), mkdir(outSubDir); end

saveS.fAll             = fAll;
saveS.cAll             = cAll;
saveS.aAll             = aAll;
saveS.bAll             = bAll;
saveS.bandIdxAll       = bandIdxAll;       % indexes into procEdgesVec
saveS.dispBandIdxAll   = dispBandIdxAll;   % indexes into dispEdgesVec
saveS.hFIR             = hFIR;
saveS.poleStats        = poleStats;        % per display band
saveS.globalStats      = globalStats;
saveS.refAI            = refAI;
saveS.nRMSE_IIR        = nRMSE_IIR;
saveS.nRMSE_final      = nRMSE_final;
saveS.dispEdgesVec     = dispEdgesVec;
saveS.procEdgesVec     = procEdgesVec;
saveS.fTransition      = fTransition;
saveS.fs               = fs;
saveS.filename         = filename;
saveS.bandMode         = bandMode;
saveS.procBandWidth_Hz = procBandWidth_Hz;
saveS.V_room           = V_room;

savePath = fullfile(outSubDir, [filename, '_polymaxRIR.mat']);
save(savePath, '-struct', 'saveS');
fprintf('MAT saved   : %s\n', savePath);

outSynNorm = outSyn / (max(abs(outSyn)) + eps);
audiowrite(fullfile(outSubDir, [filename, '_polymaxRIR_SYNTH.wav']), outSynNorm, fs);

% --- Save figures as .fig files -------------------------------------------
savefig(hFig_spectrum,  fullfile(outSubDir, [filename, '_spectrum.fig']));
savefig(hFig_bandStats, fullfile(outSubDir, [filename, '_bandStats.fig']));
savefig(hFig_LFcomp,    fullfile(outSubDir, [filename, '_LFcomparison.fig']));
fprintf('FIG saved   : %s_{spectrum,bandStats,LFcomparison}.fig\n', ...
    fullfile(outSubDir, filename));

isKept = nRMSE_final <= nRMSE_time_keep_threshold;
fprintf('\nnRMSE_final = %.4g  |  threshold = %.4g  ->  %s\n', ...
    nRMSE_final, nRMSE_time_keep_threshold, ...
    conditional(isKept, 'KEPT', 'DISCARDED'));

% --- Record in batch summary ----------------------------------------------
batchSummary.filename{iFile}        = filename;
batchSummary.nRMSE_IIR(iFile)       = nRMSE_IIR;
batchSummary.nRMSE_final(iFile)     = nRMSE_final;
batchSummary.totalPoles(iFile)      = cumPoles;
batchSummary.fTransition(iFile)     = fTransition;
batchSummary.T60global_mean(iFile)  = globalStats.T60global_mean;
batchSummary.isKept(iFile)          = isKept;
batchSummary.errorMsg{iFile}        = '';

catch ME
    % Per-file error: log and continue with the next file
    fprintf('\n  ERROR on file ''%s'' : %s\n  Skipping.\n', filename, ME.message);
    batchSummary.filename{iFile} = filename;
    batchSummary.errorMsg{iFile} = ME.message;
end   % try

close all;   % close figures before next file to avoid accumulation

end   % for iFile

%% ========================= BATCH SUMMARY =================================

fprintf('\n\n=== BATCH COMPLETE : %d/%d files processed ===\n', ...
    nnz(~cellfun(@isempty, batchSummary.filename)), NFiles);
fprintf('\n%-40s  %8s  %8s  %7s  %6s  %s\n', ...
    'Filename', 'nRMSE_IIR', 'nRMSE_fin', 'Poles', 'T60(s)', 'Status');
fprintf('%s\n', repmat('-', 1, 85));
for iFile = 1 : NFiles
    if isempty(batchSummary.errorMsg{iFile})
        status = conditional(batchSummary.isKept(iFile), 'KEPT', 'DISCARDED');
        fprintf('%-40s  %8.4g  %8.4g  %7d  %6.2f  %s\n', ...
            batchSummary.filename{iFile}, ...
            batchSummary.nRMSE_IIR(iFile), ...
            batchSummary.nRMSE_final(iFile), ...
            batchSummary.totalPoles(iFile), ...
            batchSummary.T60global_mean(iFile), ...
            status);
    else
        fprintf('%-40s  %s\n', batchSummary.filename{iFile}, ...
            ['ERROR: ', batchSummary.errorMsg{iFile}]);
    end
end

% Save batch summary alongside the individual results
batchSummary.inputFolder   = inputFolder;
batchSummary.outSubDir     = outSubDir;
batchSummary.processedDate = datestr(now);  %#ok<TNOW1,DATST>
batchSavePath = fullfile(outSubDir, 'batchSummary.mat');
save(batchSavePath, 'batchSummary');
fprintf('\nBatch summary saved : %s\n', batchSavePath);


%% ########################################################################
%%  LOCAL FUNCTIONS
%%
%%  Section order:
%%    1.  Band utilities        (buildISOBands, buildLinearBands,
%%                               selectRegime, assignToDisplayBands,
%%                               computeDisplayBandStats)
%%    2.  Parameter selection   (autoPolymaxParams, interpT60band)
%%    3.  Pole estimation       (polymaxWrapper, peakPickWrapper)
%%    4.  Amplitude solve       (leastSqWeights, rejectGainOutliers,
%%                               gainOutlierMask, globalBandedLS)
%%    5.  Spectrum / error      (spectrumBuild, mergeRejectClosePoles,
%%                               normalisedRMSE_time)
%%    6.  Statistics            (computeBandStats, emptyBandStats,
%%                               computeGlobalStats)
%%    7.  Output                (printAllStats, plotBandStats,
%%                               plotLFComparison)
%%    8.  Signal utilities      (takefft, designResidualFIR_LS)
%%    9.  Misc                  (conditional)
%% ########################################################################

% =========================================================================
%  1. Band utilities
% =========================================================================

function [edges, fc_iso] = buildISOBands(fLow, fHigh, bandsPerOctave)
%BUILDISOBANDS  Exact ISO 266 / IEC 61260 octave and third-octave band edges.
%
%   [EDGES, FC_ISO] = BUILDISOBANDS(FLOW, FHIGH, BANDSPEROCTAVE)
%
%   Returns band edges that conform strictly to ISO 266 / IEC 61260.
%   Band-edge frequencies are the exact geometric means between consecutive
%   ISO centre frequencies:
%
%       fc(n) = 1000 * 2^(n / B)     [ISO base-2 series, B = bandsPerOctave]
%       lo(n) = fc(n) / 2^(1/(2B))   [lower edge of band n]
%       hi(n) = fc(n) * 2^(1/(2B))   [upper edge of band n]
%
%   Because consecutive ISO bands are perfectly contiguous (hi(n) = lo(n+1)
%   exactly), the output EDGES is a monotonically increasing vector with no
%   gaps between bands.
%
%   KEY DESIGN DECISION
%   fLow and fHigh are used only to SELECT which ISO bands to include.
%   The band edges themselves are NEVER adjusted to match fLow / fHigh.
%   This guarantees that every band in EDGES has standard ISO width and
%   that the analysis range exactly matches integer multiples of the band
%   step — which is the correct interpretation of ISO 266.
%
%   INPUTS
%     fLow, fHigh     : desired frequency range [Hz] — used as selectors only
%     bandsPerOctave  : 1 (octave) or 3 (third-octave)
%
%   OUTPUTS
%     edges   : ISO band edge vector [Hz], length = nBands + 1
%     fc_iso  : ISO centre frequencies of the selected bands [Hz]

halfStep = 2^(1 / (2*bandsPerOctave));   % ratio from centre to band edge

% Generate a window of ISO band indices wide enough to cover [fLow, fHigh].
% We add 1 on each side so that any band partially overlapping the range
% is captured before the overlap test below filters it.
nMin = floor(bandsPerOctave * log2(fLow  / 1000)) - 1;
nMax = ceil( bandsPerOctave * log2(fHigh / 1000)) + 1;
n    = nMin : nMax;

% Exact ISO centre frequencies and their band edges
fc = 1000 * 2.^(n / bandsPerOctave);
lo = fc / halfStep;
hi = fc * halfStep;

% Retain only bands whose centre lies within [fLow, fHigh].
% Using centre-based selection (not edge-based) ensures that a band is
% included only when it is genuinely within the requested analysis range,
% even if its edges extend slightly beyond fLow / fHigh.
mask   = (fc >= fLow) & (fc <= fHigh);
lo     = lo(mask);
hi     = hi(mask);
fc_iso = fc(mask);

if isempty(lo)
    % Fallback: return a single band spanning [fLow, fHigh] if no ISO
    % centre falls within the specified range (e.g., very narrow range).
    warning('polymaxRIR:noISOBands', ...
        'No ISO band centre found in [%.1f, %.1f] Hz. Using a single band.', ...
        fLow, fHigh);
    edges  = [fLow, fHigh];
    fc_iso = sqrt(fLow * fHigh);
    return;
end

% Build the edge vector.  Because consecutive ISO bands are contiguous by
% definition (hi(k) = fc(k)*halfStep = fc(k+1)/halfStep = lo(k+1)), the
% edge vector is simply the lower edges of all bands followed by the upper
% edge of the last band.  This avoids relying on floating-point equality
% inside unique(), which can silently produce duplicate edges when hi(k) and
% lo(k+1) differ by one ULP.
edges  = [lo(:).', hi(end)];
fc_iso = fc(mask);
end

% -------------------------------------------------------------------------

function regime = selectRegime(fcBand, fTransition)
%SELECTREGIME  Return 'LF' or 'HF' based on band centre vs transition freq.

if fcBand < fTransition
    regime = 'LF';
else
    regime = 'HF';
end
end

% -------------------------------------------------------------------------

function edges = buildLinearBands(fLow, fHigh, bwHz, fv)
%BUILDLINEARBANDS  Narrow linear processing bands of fixed width BWHZ.
%
%   Produces edges [fLow, fLow+bwHz, fLow+2*bwHz, ...] clipped so that
%   the last edge does not exceed fHigh and all edges remain within the
%   usable FFT range [fv(2), fv(end)].  A partial band at the top is
%   always included so that no part of [fLow, fHigh] is left uncovered.

edges = fLow : bwHz : fHigh;
if edges(end) < fHigh
    edges(end+1) = fHigh;   % add partial band at top if needed
end
% Clip to FFT range
edges(1)   = max(edges(1),   fv(2));
edges(end) = min(edges(end), fv(end));
edges      = unique(edges);
end

% -------------------------------------------------------------------------

function dispIdx = assignToDisplayBands(fAll, dispEdgesVec)
%ASSIGNTODISPLAYBANDS  Map each pole to its ISO display band by frequency.
%
%   Returns DISPIDX (same length as FALL) where DISPIDX(p) = k means
%   pole p falls in display band k, i.e. dispEdgesVec(k) <= fAll(p) < dispEdgesVec(k+1).
%   Poles outside the display range receive index 0.

NDisp  = numel(dispEdgesVec) - 1;
dispIdx = zeros(numel(fAll), 1, 'int32');
for k = 1 : NDisp
    mask = (fAll >= dispEdgesVec(k)) & (fAll < dispEdgesVec(k+1));
    dispIdx(mask) = k;
end
% Include any poles sitting exactly on the last edge
dispIdx(fAll == dispEdgesVec(end) & dispIdx == 0) = NDisp;
end

% -------------------------------------------------------------------------

function poleStats = computeDisplayBandStats(fAll, cAll, dispEdgesVec, ...
    dispCentresISO, dispBandIdxAll, fTransition)
%COMPUTEDISPLAYBANDSTATS  Build per-display-band statistics from all poles.
%
%   Groups poles by their display band index, computes statistics for each
%   band, and sets regime ('LF'/'HF') based on whether the band centre is
%   below or above fTransition.  Preserves the per-pole T60 spread (stdT60)
%   correctly since cAll at this point still contains the pre-homogenisation
%   values stored per processing band (all poles in a processing band share
%   the same c, but processing bands are narrower than display bands so the
%   spread across processing bands within one display band is meaningful).

NDisp = numel(dispEdgesVec) - 1;

poleStats_template = struct( ...
    'flow', NaN, 'fhigh', NaN, 'fcBand', NaN, 'N', 0, ...
    'densityRatio', NaN, 'modalDensity_Hz', NaN, ...
    'meanF', NaN, 'stdF', NaN, ...
    'meanSpacing_Hz', NaN, 'stdSpacing_Hz', NaN, 'cvSpacing', NaN, ...
    'meanT60', NaN, 'stdT60', NaN, 'cvT60', NaN, ...
    'meanC', NaN, 'stdC', NaN, 'regime', 'LF');
poleStats = repmat(poleStats_template, 1, NDisp);

for k = 1 : NDisp
    fc     = dispCentresISO(k);
    regime = selectRegime(fc, fTransition);
    idx    = (dispBandIdxAll == k);

    if ~any(idx)
        poleStats(k) = emptyBandStats(dispEdgesVec(k), dispEdgesVec(k+1), fc, regime);
        continue;
    end

    poleStats(k) = computeBandStats( ...
        fAll(idx), cAll(idx), dispEdgesVec(k), dispEdgesVec(k+1), fc, regime);
end
end
% =========================================================================

function [NordMin, NordMax, ordStep, T60cut, Nexpected] = ...
    autoPolymaxParams(fcBand, bwBand, Nbins, V_room, c_sound, T60est, ...
    maxOrderLF, maxIters, T60cutMult, T60cutFloor)
%AUTOPOLYMAXPARAMS  Automatically select PolyMax order for one LF band.
%
%   Uses Weyl's modal density formula to estimate the expected number of
%   modes in the band and sizes the polynomial order accordingly:
%
%     N_expected = round( (4π V fc²/c³) × BW )
%     NordMax    = min( 6 × N_expected,  Nbins/4,  maxOrderLF )
%     T60cut     = max( T60cutFloor,  T60cutMult × T60est )
%
%   INPUTS
%     fcBand, bwBand : band centre [Hz] and bandwidth [Hz]
%     Nbins          : number of FFT bins in this band
%     V_room, c_sound: room volume [m³] and speed of sound [m/s]
%     T60est         : estimated T60 for this band [s]
%     maxOrderLF     : hard ceiling on NordMax
%     maxIters       : target number of order steps (sets ordStep)
%     T60cutMult     : T60cut multiplier (e.g. 3.0)
%     T60cutFloor    : minimum T60cut [s]
%
%   OUTPUTS
%     NordMin, NordMax, ordStep : order sweep parameters
%     T60cut                    : T60 rejection ceiling [s]
%     Nexpected                 : estimated mode count in band

% Expected mode count from Weyl's approximation: dn/df = 4π V fc² / c³
dndF      = 4*pi * V_room * fcBand^2 / c_sound^3;
Nexpected = max(1, round(dndF * bwBand));

NordMin      = 2;
NordMax_bins = floor(Nbins / 4);

% Cap T60est to prevent runaway T60cut on pathological inputs
T60cut = max(T60cutFloor, T60cutMult * min(T60est, 120));

NordMax = min([maxOrderLF, round(6.0 * Nexpected), NordMax_bins]);
NordMax = min(NordMax, 900);   % hard cap: keeps companion matrix tractable
NordMax = max(NordMax, NordMin + 20);

ordStep = max(1, round((NordMax - NordMin) / maxIters));
end

% -------------------------------------------------------------------------

function T60est = interpT60band(fcBand, fcBandsAI, T20bandsAI, cFallback)
%INTERPT60BAND  Interpolate T60 estimate at an arbitrary band centre.
%
%   Interpolates log-linearly through the per-band T30 values from
%   acoustic_indices, extrapolating at the edges.  Falls back to
%   3*ln(10)/cFallback if fewer than two valid T30 values are available.

T60all = T20bandsAI(:);     % T20 is an ISO estimate of T60; no factor of 2
idxV   = ~isnan(T60all) & (T60all > 0);

if nnz(idxV) >= 2
    T60est = interp1(log(fcBandsAI(idxV)), T60all(idxV), ...
        log(max(fcBand, 1)), 'linear', 'extrap');
    T60est = max(0.1, T60est);
elseif nnz(idxV) == 1
    T60est = T60all(find(idxV, 1));
else
    T60est = 3*log(10) / cFallback;
end
end

% =========================================================================
%  3. Pole estimation
% =========================================================================

function [fBand, cBand] = polymaxWrapper(fftBand, fvBand, ...
    NordMin, NordMax, ordStep, tolf, told, T60cut, Nexpected, doLivePlot)
%POLYMAXWRAPPER  Run POLYMAX for one LF band with automatic order retry.
%
%   Calls polymax.m with NoCluster = 1 (stability-diagram poles returned
%   directly) and retries with a progressively higher NordMax if the result
%   contains fewer than 30 % of the expected mode count.  Up to three
%   attempts are made before accepting the best result found.
%
%   When DOLIVEPLOT is true, a figure is opened before each call so that
%   the stability diagram renders in a clean single-axes window.

Nexpected_min = max(1, floor(0.3 * Nexpected));
maxRetries    = 3;
NordMax_cur   = NordMax;
NordMax_hard  = min(floor(numel(fvBand) / 4), 2000);
fBand = [];
cBand = [];

for iTry = 1 : maxRetries

    if doLivePlot
        figure('Name', sprintf('PolyMax [%.1f – %.1f Hz]', fvBand(1), fvBand(end)));
    end

    try
        [fTry, cTry, ~] = polymax(fftBand(:), fvBand(:), ...
            NordMin, NordMax_cur, ordStep, tolf, told, ...
            0, T60cut, 0, 1, double(doLivePlot));
    catch ME
        warning('polymaxRIR:polymaxFailed', '%s', ME.message);
        fTry = [];
        cTry = [];
    end

    if numel(fTry) >= Nexpected_min
        fBand = fTry;
        cBand = cTry;
        return;
    end

    % Keep the best result seen so far
    if numel(fTry) > numel(fBand)
        fBand = fTry;
        cBand = cTry;
    end

    % Retry with a higher order, up to the hard bin-count ceiling
    if iTry < maxRetries && NordMax_cur < NordMax_hard
        NordMax_cur = min(NordMax_hard, NordMax_cur + 20);
    else
        return;
    end
end
end

% -------------------------------------------------------------------------

function [fPeaks, cPeaks] = peakPickWrapper(fftBand, fvBand, df, ...
    minProm_dB, minDist_Hz, cFixed)
%PEAKPICKWRAPPER  HF peak detection with per-band Schroeder damping.
%
%   Detects peaks in |H(f)| using findpeaks.  Damping is set uniformly to
%   CFIXED = 3*ln(10)/T60est from the Schroeder T30 — the only reliable
%   damping estimate above fTransition where modal overlap invalidates any
%   per-peak bandwidth method.
%
%   Prominence is evaluated relative to the local spectral floor (median
%   level), converted to linear scale before calling findpeaks.

fftBand = fftBand(:);
fvBand  = fvBand(:);
magBand = abs(fftBand);
Nb      = numel(fvBand);

if Nb < 3 || all(magBand == 0)
    fPeaks = [];
    cPeaks = [];
    return;
end

dB_band      = 20*log10(magBand + eps);
minDistSamp  = max(1, round(minDist_Hz / max(df, eps)));
localFloor   = median(dB_band);                     % [dB] spectral floor estimate
floorLin     = 10^(localFloor / 20);
minProm_lin  = max(eps, 10^((localFloor + minProm_dB) / 20) - floorLin);

try
    [~, locs] = findpeaks(magBand, ...
        'MinPeakProminence', minProm_lin, ...
        'MinPeakDistance',   minDistSamp);
catch
    % Fallback: simple discrete local-maxima without prominence
    locs = find((magBand(2:end-1) > magBand(1:end-2)) & ...
                (magBand(2:end-1) >= magBand(3:end))) + 1;
end

if isempty(locs)
    fPeaks = [];
    cPeaks = [];
    return;
end

fPeaks = fvBand(locs);
cPeaks = cFixed * ones(size(fPeaks));
end

% =========================================================================
%  4. Amplitude solve
% =========================================================================

function [aBand, bBand] = leastSqWeights(omvBand, Npoles, fBand, cBand, fftRef)
%LEASTSQWEIGHTS  Least-squares residue amplitudes for a modal partial-fraction model.
%
%   The frequency-domain model is the sum of conjugate pole pairs:
%
%     H(jω) = Σ_p  [ r_p / (jω - λ_p)  +  r_p* / (jω - λ_p*) ]
%
%   where λ_p = -c_p + j√(ω_p² - c_p²) and the complex residue is
%   r_p = a_p + j·b_p.  Substituting and separating real parts yields
%   a linear system in [a_p; b_p] that is solved here.
%
%   Tikhonov regularisation is amplitude-invariant: the parameter scales
%   with the mean diagonal of R'R so results are independent of signal level.

Lfv   = numel(omvBand);
fBand = fBand(:).';
cBand = cBand(:).';

rad   = max((2*pi*fBand).^2 - cBand.^2, 1e-12);
poles = -cBand + 1j*sqrt(rad);       % 1×Npoles

% Build basis matrices Ra, Rb via broadcasting (Lfv×Npoles each).
% jomv is Lfv×1, poles is 1×Npoles → D is Lfv×Npoles without a loop.
jomv = 1j * omvBand(:);              % Lfv×1
D    = jomv - poles;                 % Lfv×Npoles:  jω - λ_p
Dc   = jomv - conj(poles);          % Lfv×Npoles:  jω - λ_p*
Ra   = 1./D  + 1./Dc;               % Lfv×Npoles
Rb   = 1j./D - 1j./Dc;             % Lfv×Npoles

R   = [Ra, Rb];
RtR = R' * R;

% Amplitude-invariant regularisation
eps_rel = 1e-8;
reg     = eps_rel * trace(RtR) / (2*Npoles) * eye(2*Npoles);

solb = (RtR + reg) \ (R' * fftRef(:));
aT   = solb(1:Npoles);
bT   = solb(Npoles + 1 : end);

% Extract physical (real-valued) amplitude components
aBand = real(aT) - imag(bT);
bBand = real(bT) + imag(aT);
end

% -------------------------------------------------------------------------

function [aOut, bOut, fOut, cOut] = rejectGainOutliers( ...
    aIn, bIn, fIn, cIn, omvBand, fftBand, method, zthresh, maxIters, minModes)
%REJECTGAINOUTLIERS  Iteratively remove poles with anomalous residue amplitudes.
%
%   Flags poles whose real or imaginary amplitude deviates more than ZTHRESH
%   robust Z-scores from the band median, removes them, and re-solves for
%   amplitudes.  Iterates up to MAXITERS times or until no outliers remain.
%   Never reduces the pole count below MINMODES.

aOut = aIn(:);
bOut = bIn(:);
fOut = fIn(:);
cOut = cIn(:);

for it = 1 : maxIters
    if numel(aOut) <= minModes, break; end
    mask = gainOutlierMask(aOut, bOut, method, zthresh);
    if ~any(mask) || nnz(~mask) < minModes, break; end
    aOut = aOut(~mask);
    bOut = bOut(~mask);
    fOut = fOut(~mask);
    cOut = cOut(~mask);
    [aOut, bOut] = leastSqWeights(omvBand, numel(fOut), fOut, cOut, fftBand);
end
end

% -------------------------------------------------------------------------

function mask = gainOutlierMask(a, b, method, zthresh)
%GAINOUTLIERMASK  Boolean mask: true where residues are outliers.
%
%   method = 'mad':  robust Z-score using median ± 1.4826 × MAD.
%   method = 'std':  classical Z-score using mean ± std.

a = a(:);
b = b(:);

switch lower(method)
    case 'mad'
        ma = median(a);  mb = median(b);
        sa = max(1.4826 * median(abs(a - ma)), eps);
        sb = max(1.4826 * median(abs(b - mb)), eps);
    case 'std'
        ma = mean(a);  mb = mean(b);
        sa = std(a) + eps;
        sb = std(b) + eps;
    otherwise
        error('polymaxRIR:unknownMethod', ...
            'Unknown gainOutlierMethod ''%s''.  Use ''mad'' or ''std''.', method);
end

mask = (abs((a - ma) ./ sa) > zthresh) | (abs((b - mb) ./ sb) > zthresh);
end

% -------------------------------------------------------------------------

function [aNew, bNew] = globalBandedLS( ...
    spcTarget, fv, fAll, cAll, bandIdxAll, bandEdgesVec, overlapFrac)
%GLOBALBANDEDLS  Global amplitude re-solve with per-band overlapping windows.
%
%   Re-fits residue amplitudes for all poles simultaneously, processing one
%   band at a time but extending each LS window by OVERLAPFRAC on each side
%   to reduce edge artefacts at band boundaries.

NBands = numel(bandEdgesVec) - 1;
aNew   = zeros(size(fAll));
bNew   = zeros(size(fAll));

for ib = 1 : NBands
    flow  = bandEdgesVec(ib);
    fhigh = bandEdgesVec(ib + 1);
    bw    = fhigh - flow;

    idxP = (bandIdxAll == ib);
    Np   = nnz(idxP);
    if Np == 0, continue; end

    [~, is] = min(abs(fv - max(fv(1),   flow  - overlapFrac*bw)));
    [~, ie] = min(abs(fv - min(fv(end), fhigh + overlapFrac*bw)));
    ie      = max(ie, is + 1);

    omvSolve = 2*pi*fv(is:ie);
    fftSolve = spcTarget(is:ie);

    [aBand, bBand] = leastSqWeights(omvSolve(:), Np, ...
        fAll(idxP), cAll(idxP), fftSolve(:));

    idx        = find(idxP);
    aNew(idx)  = aBand;
    bNew(idx)  = bBand;
end
end

% =========================================================================
%  5. Spectrum / error
% =========================================================================

function h = synthesizeIR(aVec, bVec, cVec, fVec, fs, Nout)
%SYNTHESIZEIR  Direct time-domain modal synthesis — no IFFT, no wrap-around.
%
%   h(n) = 2·Re[ Σ_p r_p · exp(λ_p · n/fs) ]
%   r_p  = a_p + j·b_p
%   λ_p  = -c_p + j·√( (2π f_p)² - c_p² )
%
%   An IFFT-based reconstruction treats the spectrum as periodic with
%   period T = N/fs.  Any pole whose T60 exceeds T has its tail wrapped
%   back onto the end of the block, producing a spurious burst.
%   This function evaluates the modal sum directly at each time sample
%   so every mode decays correctly to zero with no periodicity assumption.
%   Poles are processed in chunks to keep memory usage bounded.

aVec = aVec(:);  bVec = bVec(:);  cVec = cVec(:);  fVec = fVec(:);
Np     = numel(aVec);
rad    = max((2*pi*fVec).^2 - cVec.^2, 1e-12);
lambda = -cVec + 1j*sqrt(rad);   % complex pole positions [rad/s]
resids = aVec + 1j*bVec;

n = (0 : Nout-1);   % 1×Nout time-sample index

% Chunk size: keep each partial matrix ≤ ~50 MB (complex128 = 16 bytes)
chunkSize = max(1, floor(50e6 / (Nout * 16)));

h = zeros(Nout, 1);
for p = 1 : chunkSize : Np
    idx = p : min(p + chunkSize - 1, Np);
    h = h + 2 * real(sum(resids(idx) .* exp(lambda(idx) .* n / fs), 1)).';
end
end

% -------------------------------------------------------------------------

function spc = spectrumBuildSlice(aVec, bVec, cVec, fVec, omv)
%SPECTRUMBUILDSLICE  Modal spectrum for one narrow band (small Np, no chunking).
%
%   Identical algebra to spectrumBuild but operates on a short omv slice
%   with only the Np poles belonging to one processing band, so the
%   Np×Nbins_slice matrix is always small and memory-safe.  Used inside
%   the iterative refinement loop to accumulate spcSyn_cur band-by-band.

aVec = aVec(:);  bVec = bVec(:);  cVec = cVec(:);  fVec = fVec(:);
omv  = omv(:).';
rad    = max((2*pi*fVec).^2 - cVec.^2, 1e-12);
poles  = -cVec + 1j*sqrt(rad);
resids = aVec  + 1j*bVec;
D   = 1j*omv - poles;
Dc  = 1j*omv - conj(poles);
spc = sum(resids./D + conj(resids)./Dc, 1);
spc = spc(:);
end

% -------------------------------------------------------------------------

function spc = spectrumBuild(aVec, bVec, cVec, fVec, omv)
%SPECTRUMBUILD  Evaluate the modal partial-fraction spectrum at OMV.
%
%   Reconstructs H(jω) = Σ [ r_p/(jω−λ_p) + r_p*/(jω−λ_p*) ] at each
%   frequency in OMV, where λ_p = −c_p + j·√(ωp²−cp²) and r_p = a_p + j·b_p.
%
%   Implementation: poles are processed in chunks so that each broadcast
%   (Nchunk×Nfft) matrix fits comfortably in memory.  This is far faster than
%   a per-pole loop but avoids allocating a full Npoles×Nfft matrix which
%   would exhaust RAM for large pole counts.

aVec = aVec(:);
bVec = bVec(:);
cVec = cVec(:);
fVec = fVec(:);
omv  = omv(:).';           % 1×Nfft

Np     = numel(aVec);
Nfft   = numel(omv);
rad    = max((2*pi*fVec).^2 - cVec.^2, 1e-12);
poles  = -cVec + 1j*sqrt(rad);   % Np×1
resids = aVec + 1j*bVec;         % Np×1

% Chunk size: target ~50 MB per chunk (complex double = 16 bytes)
chunkSize = max(1, floor(50e6 / (Nfft * 16)));

spc = zeros(1, Nfft);
for p = 1 : chunkSize : Np
    idx = p : min(p + chunkSize - 1, Np);
    r  = resids(idx);      % chunk×1
    pl = poles(idx);       % chunk×1
    % Broadcast: (chunk×1) ./ (1×Nfft - chunk×1) → chunk×Nfft
    D  = 1j*omv - pl;      % chunk×Nfft
    Dc = 1j*omv - conj(pl);
    spc = spc + sum(r./D + conj(r)./Dc, 1);
end
spc = spc(:);
end

% -------------------------------------------------------------------------

function [fKeep, cKeep] = mergeRejectClosePoles(fNew, cNew, fExist, tolHz)
%MERGEREJECTCLOSEPOLES  Discard new poles that are within TOLHZ of existing ones.
%
%   Used during residual refinement to prevent a mode from being represented
%   by two nearly identical poles across successive iterations.

fNew   = fNew(:);
cNew   = cNew(:);
fExist = fExist(:);

keep = true(size(fNew));
for k = 1 : numel(fNew)
    if ~isempty(fExist) && any(abs(fExist - fNew(k)) <= tolHz)
        keep(k) = false;
    end
end

fKeep = fNew(keep);
cKeep = cNew(keep);
end

% -------------------------------------------------------------------------

function nrmse = normalisedRMSE_time(xRef, xSyn)
%NORMALISEDRMSE_TIME  Time-domain normalised RMSE.
%
%   nRMSE = RMS(xSyn - xRef) / max(|xRef|)
%   Computed over the shorter of the two signals.

xRef  = xRef(:);
xSyn  = xSyn(:);
L     = min(numel(xRef), numel(xSyn));
nrmse = sqrt(mean((xSyn(1:L) - xRef(1:L)).^2)) / (max(abs(xRef(1:L))) + eps);
end

% =========================================================================
%  6. Statistics
% =========================================================================

function stats = computeBandStats(fBand, cBand, flow, fhigh, fcBand, regime)
%COMPUTEBANDSTATS  Per-band modal statistics struct.
%
%   Computes mode count, modal density, mean/std of inter-mode spacing, T60,
%   and damping for the poles in FBAND/CBAND.

fBand = fBand(:);
cBand = cBand(:);
bw    = fhigh - flow;
N     = numel(fBand);

stats.flow    = flow;
stats.fhigh   = fhigh;
stats.fcBand  = fcBand;
stats.N       = N;
stats.regime  = regime;
stats.densityRatio = NaN;   % filled in by computeGlobalStats

if N == 0
    [stats.modalDensity_Hz, stats.meanF, stats.stdF, ...
     stats.meanSpacing_Hz,  stats.stdSpacing_Hz, stats.cvSpacing, ...
     stats.meanT60, stats.stdT60, stats.cvT60, ...
     stats.meanC,   stats.stdC] = deal(NaN);
    return;
end

stats.modalDensity_Hz = N / bw;
stats.meanF = mean(fBand);
stats.stdF  = std(fBand);

if N >= 2
    sp = diff(sort(fBand));
    stats.meanSpacing_Hz = mean(sp);
    stats.stdSpacing_Hz  = std(sp);
    stats.cvSpacing      = stats.stdSpacing_Hz / (stats.meanSpacing_Hz + eps);
else
    stats.meanSpacing_Hz = bw;
    stats.stdSpacing_Hz  = 0;
    stats.cvSpacing      = 0;
end

T60vec     = 3*log(10) ./ max(cBand, 1e-10);
stats.meanT60 = mean(T60vec);
stats.stdT60  = std(T60vec);
stats.cvT60   = stats.stdT60 / (stats.meanT60 + eps);
stats.meanC   = mean(cBand);
stats.stdC    = std(cBand);
end

% -------------------------------------------------------------------------

function stats = emptyBandStats(flow, fhigh, fcBand, regime)
%EMPTYBANDSTATS  Return a zero-pole statistics struct for a skipped band.

stats = computeBandStats([], [], flow, fhigh, fcBand, regime);
end

% -------------------------------------------------------------------------

function gStats = computeGlobalStats(poleStats, fAll, cAll, bandEdgesVec, V_room, c_sound)
%COMPUTEGLOBALSTATS  Aggregate per-band statistics into a global summary.
%
%   Computes the global mean T60, estimated Schroeder frequency, and the
%   theoretical / measured mode count below the Schroeder frequency.
%   Also annotates each poleStats entry with its density ratio relative to
%   the Weyl approximation.

NBands = numel(poleStats);
cPos   = cAll(cAll > 0);
T60global = mean(3*log(10) ./ cPos, 'omitnan');
fSch      = 2000 * sqrt(max(T60global, 0.1) / V_room);

gStats.totalPoles       = numel(fAll);
gStats.T60global_mean   = T60global;
gStats.fSchroeder       = fSch;
gStats.V_room           = V_room;
gStats.Ntheoretical_Sch = round((4*pi/3) * V_room * (fSch / c_sound)^3);
gStats.Nactual_Sch      = nnz(fAll <= fSch);

fcVec    = zeros(NBands, 1);
NvecB    = zeros(NBands, 1);
densVec  = zeros(NBands, 1);
densTheo = zeros(NBands, 1);
denRatio = zeros(NBands, 1);
mT60Vec  = zeros(NBands, 1);
sT60Vec  = zeros(NBands, 1);
cvT60Vec = zeros(NBands, 1);
mSpVec   = zeros(NBands, 1);
sSpVec   = zeros(NBands, 1);
cvSpVec  = zeros(NBands, 1);

for ib = 1 : NBands
    s    = poleStats(ib);
    fc   = s.fcBand;
    dndF = 4*pi * V_room * fc^2 / c_sound^3;   % Weyl modal density [modes/Hz]

    fcVec(ib)    = fc;
    NvecB(ib)    = s.N;
    densVec(ib)  = s.modalDensity_Hz;
    densTheo(ib) = dndF;
    if dndF > 0
        denRatio(ib) = s.modalDensity_Hz / dndF;
    end
    mT60Vec(ib)  = s.meanT60;
    sT60Vec(ib)  = s.stdT60;
    cvT60Vec(ib) = s.cvT60;
    mSpVec(ib)   = s.meanSpacing_Hz;
    sSpVec(ib)   = s.stdSpacing_Hz;
    cvSpVec(ib)  = s.cvSpacing;
    poleStats(ib).densityRatio = denRatio(ib);  %#ok<NASGU>  (caller sees poleStats, not ib-indexed)
end

gStats.fcVec    = fcVec;
gStats.NvecB    = NvecB;
gStats.densVec  = densVec;
gStats.densTheo = densTheo;
gStats.denRatio = denRatio;
gStats.mT60Vec  = mT60Vec;
gStats.sT60Vec  = sT60Vec;
gStats.cvT60Vec = cvT60Vec;
gStats.mSpVec   = mSpVec;
gStats.sSpVec   = sSpVec;
gStats.cvSpVec  = cvSpVec;
end

% =========================================================================
%  7. Output
% =========================================================================

function printAllStats(gStats, poleStats, verbose)
%PRINTALLSTATS  Print global and per-band statistics to the console.

fprintf('\n--- GLOBAL SUMMARY ---\n');
fprintf('  Total poles              : %d\n',      gStats.totalPoles);
fprintf('  Mean T60 (all poles)     : %.3f s\n',  gStats.T60global_mean);
fprintf('  Schroeder frequency      : %.1f Hz\n', gStats.fSchroeder);
fprintf('  Modes <= f_Sch (actual)  : %d\n',      gStats.Nactual_Sch);
fprintf('  Modes <= f_Sch (theory)  : %d\n',      gStats.Ntheoretical_Sch);

if ~verbose, return; end

fprintf('\n%4s %3s %8s %8s %5s %9s %9s %8s %8s %8s\n', ...
    'Band', 'Reg', 'flow', 'fhigh', 'N', 'dens_Hz', 'theo_Hz', ...
    'mT60', 'cvT60', 'mSpc');

for ib = 1 : numel(poleStats)
    s  = poleStats(ib);
    th = gStats.densTheo(ib);
    fprintf('%4d %3s %8.1f %8.1f %5d %9.4f %9.4f %8.3f %8.3f %8.3f\n', ...
        ib, s.regime, s.flow, s.fhigh, s.N, s.modalDensity_Hz, th, ...
        s.meanT60, s.cvT60, s.meanSpacing_Hz);
end
end

% -------------------------------------------------------------------------

function plotBandStats(poleStats, gStats, NBands, fTransition)
%PLOTBANDSTATS  Three-panel per-band modal statistics figure.
%
%   Panels:
%     1 — Mode count per band (blue = LF, orange = HF)
%     2 — Mean ± 1σ T60 per band (log-frequency axis)
%     3 — Mean ± 1σ inter-mode spacing (log-frequency axis)

fcVec = gStats.fcVec;
NvecB = gStats.NvecB;
mT60  = gStats.mT60Vec;
sT60  = gStats.sT60Vec;
mSp   = gStats.mSpVec;
sSp   = gStats.sSpVec;
isHF  = fcVec >= fTransition;
xLab  = arrayfun(@(f) sprintf('%.0f', f), fcVec, 'UniformOutput', false);

subplot(1, 3, 1);
b = bar(1:NBands, NvecB, 'EdgeColor', 'none');
b.FaceColor = 'flat';
for ib = 1 : NBands
    b.CData(ib, :) = conditional(isHF(ib), [0.85 0.45 0.20], [0.25 0.45 0.75]);
end
set(gca, 'XTick', 1:NBands, 'XTickLabel', xLab);
xtickangle(60);
ylabel('Mode count');
title('Modes per band  (blue = LF, orange = HF)');

subplot(1, 3, 2);
errorbar(fcVec, mT60, sT60, 'o-', 'Color', [0.8 0.2 0.15], ...
    'LineWidth', 1.2, 'MarkerSize', 5);
xline(fTransition, 'g--', 'f_{Sch}', 'LineWidth', 1);
set(gca, 'XScale', 'log');
xlabel('f_c (Hz)');  ylabel('T_{60} (s)');
title('Mean \pm 1\sigma T_{60} per band');

subplot(1, 3, 3);
errorbar(fcVec, mSp, sSp, 's-', 'Color', [0.35 0.25 0.65], ...
    'LineWidth', 1.2, 'MarkerSize', 5);
xline(fTransition, 'g--', 'f_{Sch}', 'LineWidth', 1);
set(gca, 'XScale', 'log');
xlabel('f_c (Hz)');  ylabel('\Deltaf (Hz)');
title('Mean \pm 1\sigma inter-mode spacing');

sgtitle('Per-band pole statistics — polymaxRIR');
end

% -------------------------------------------------------------------------

function plotLFComparison(poleStats, refAI, bandEdgesVec, bandCentresISO, fTransition)
%PLOTLFCOMPARISON  Compare PolyMax T60 against acoustic_indices T30 for LF bands.
%
%   Produces a single-panel plot showing, for every analysis band below
%   fTransition (the Schroeder frequency):
%
%     • PolyMax T60  — mean ± 1σ computed from the individual pole damping
%                      coefficients in each band:  T60_p = 3·ln10 / c_p.
%                      Error bars span [mean − std, mean + std].
%                      Plotted as filled circles with vertical error bars.
%
%     • Acoustic-indices T30  — the reverberation time estimated from the
%                      Schroeder integral by ISO 3382 (T30 fitted over -5 to
%                      -35 dB, extrapolated to 60 dB).  T30 is itself an
%                      estimate of T60 and is plotted directly as such.
%
%   The two traces use the same ISO centre frequencies, making direct
%   comparison straightforward.  Bands with no PolyMax poles (N = 0) or
%   NaN acoustic-indices T30 are still plotted if the other source has data,
%   using open markers to signal the missing estimate.
%
%   INPUTS
%     poleStats      : per-band statistics struct array from polymaxRIR
%     refAI          : results struct from acoustic_indices
%     bandEdgesVec   : ISO band edge vector [Hz]
%     bandCentresISO : ISO centre frequency of each band [Hz]
%     fTransition    : Schroeder frequency [Hz]

% --- Identify LF bands ----------------------------------------------------
NBands = numel(poleStats);
isLF   = false(1, NBands);
for ib = 1 : NBands
    if ~isnan(poleStats(ib).fcBand)
        isLF(ib) = (poleStats(ib).fcBand < fTransition);
    end
end

if ~any(isLF)
    title('No LF bands found below f_{Sch}');
    return;
end

lfIdx = find(isLF);
NLF   = numel(lfIdx);

% --- Extract PolyMax T60 statistics for LF bands -------------------------
fc_lf      = bandCentresISO(lfIdx);
pm_T60mean = nan(NLF, 1);
pm_T60std  = nan(NLF, 1);
pm_N       = zeros(NLF, 1);

for k = 1 : NLF
    ib = lfIdx(k);
    if poleStats(ib).N > 0 && ~isnan(poleStats(ib).meanT60)
        pm_T60mean(k) = poleStats(ib).meanT60;
        pm_T60std(k)  = poleStats(ib).stdT60;
        pm_N(k)       = poleStats(ib).N;
    end
end

% --- Match acoustic_indices T30 to LF band centres -----------------------
% acoustic_indices bands are at ISO centres with the same bandsPerOctave,
% but may span a different frequency range.  We find the closest AI band
% centre to each LF analysis band centre using log-frequency distance.

ai_fc  = refAI.bandCentersHz(:);
ai_T20 = refAI.T20(:);

ai_T60atLF = nan(NLF, 1);

for k = 1 : NLF
    [~, idx] = min(abs(log2(ai_fc / max(fc_lf(k), eps))));
    if ~isempty(idx) && ~isnan(ai_T20(idx)) && ai_T20(idx) > 0
        ai_T60atLF(k) = ai_T20(idx);   % T20 is an ISO 3382 estimate of T60
    end
end

% --- X-axis positions (integer indices, labels = ISO freq) ---------------
xPos = 1 : NLF;
xLab = arrayfun(@(f) sprintf('%.0f', f), fc_lf, 'UniformOutput', false);

% --- Plot -----------------------------------------------------------------
hold on;

% PolyMax T60: mean ± 1σ error bars
hasPM = ~isnan(pm_T60mean);
if any(hasPM)
    errorbar(xPos(hasPM), pm_T60mean(hasPM), pm_T60std(hasPM), ...
        'LineStyle',       'none', ...
        'Marker',          'o', ...
        'Color',           [0.18 0.40 0.75], ...
        'MarkerFaceColor', [0.18 0.40 0.75], ...
        'MarkerSize',      7, ...
        'LineWidth',       1.5, ...
        'CapSize',         8, ...
        'DisplayName',     'PolyMax  T_{60}  (mean \pm 1\sigma)');
end

% Open markers for LF bands where PolyMax found 0 poles
hasNoPM = isnan(pm_T60mean) & ~isnan(ai_T60atLF);
if any(hasNoPM)
    plot(xPos(hasNoPM), nan(nnz(hasNoPM), 1), ...
        'o', 'MarkerEdgeColor', [0.18 0.40 0.75], ...
        'MarkerFaceColor', 'none', 'MarkerSize', 7, ...
        'DisplayName', 'PolyMax  (no poles)');
end

% Acoustic-indices T30 (= T60 estimate per ISO 3382)
hasAI = ~isnan(ai_T60atLF);
if any(hasAI)
    plot(xPos(hasAI), ai_T60atLF(hasAI), ...
        's', ...
        'Color',           [0.82 0.25 0.15], ...
        'MarkerFaceColor', [0.82 0.25 0.15], ...
        'MarkerSize',      8, ...
        'LineStyle',       '--', ...
        'LineWidth',       1.1, ...
        'DisplayName',     'Acoustic indices  T_{30}  (ISO 3382)');

    % Connect with dashed line for readability
    plot(xPos(hasAI), ai_T60atLF(hasAI), '--', ...
        'Color', [0.82 0.25 0.15], 'LineWidth', 1.1, ...
        'HandleVisibility', 'off');
end

% Annotate mode count above each PolyMax marker
for k = 1 : NLF
    if pm_N(k) > 0
        text(xPos(k), pm_T60mean(k) + pm_T60std(k) + 0.03, ...
            sprintf('%d', pm_N(k)), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment',   'bottom', ...
            'FontSize', 7, 'Color', [0.18 0.40 0.75]);
    end
end

% Formatting
set(gca, 'XTick', xPos, 'XTickLabel', xLab);
xtickangle(45);
xlabel('ISO band centre frequency (Hz)');
ylabel('T_{60} (s)');
legend('show', 'Location', 'best');
grid on;
title(sprintf(['LF bands (f_c < %.0f Hz) — PolyMax T_{60} vs acoustic-indices T_{30}  (ISO 3382)\n' ...
    '(numbers above error bars = pole count per band)'], fTransition));
end

% =========================================================================
%  8. Signal utilities
% =========================================================================

function [fftinput, fv, omv] = takefft(input, fs, fftlength)
%TAKEFFT  Single-sided FFT normalised by the sampling period.
%
%   Returns the positive-frequency half of the FFT scaled by 1/fs so that
%   the result approximates the continuous Fourier transform when multiplied
%   by the frequency resolution df.

N        = max(fftlength, numel(input));
fftinput = fft(input, N) / fs;
fftinput = fftinput(1 : floor(N/2));
fv       = (0 : N-1).' * fs / N;
fv       = fv(1 : floor(N/2));
omv      = 2*pi*fv;
end

% -------------------------------------------------------------------------

function [hFIR, outFIR] = designResidualFIR_LS(targetIR, iirIR, Nfir, lambdaFIR, Nfit)
%DESIGNRESIDUALFIR_LS  Least-squares FIR residual correction filter.
%
%   Designs a length-NFIR FIR filter h such that conv(iirIR, h) ≈ targetIR
%   over the first NFIT samples.  The solve is Tikhonov-regularised with
%   parameter LAMBDAFIR to control filter energy.
%
%   Used as an optional post-processing step to compensate for systematic
%   spectral residuals that the modal IIR model cannot represent.

targetIR = targetIR(:);
iirIR    = iirIR(:);
Nfit     = min([Nfit, numel(targetIR), numel(iirIR)]);

d = targetIR(1:Nfit);
y = iirIR(1:Nfit);

% Build Toeplitz convolution matrix
Y = zeros(Nfit, Nfir);
for j = 1 : Nfir
    idxR       = j : Nfit;
    Y(idxR, j) = y(1 : numel(idxR));
end

hFIR   = (Y.'*Y + lambdaFIR*eye(Nfir)) \ (Y.'*d);
outFIR = conv(iirIR, hFIR);
end

% =========================================================================
%  9. Misc
% =========================================================================

function s = conditional(cond, a, b)
%CONDITIONAL  Inline if-else returning A when COND is true, B otherwise.
%
%   Equivalent to the ternary operator (? :) in C-family languages.
%   Avoids multi-line if/else blocks for simple scalar or string choices.

if cond
    s = a;
else
    s = b;
end
end