%==========================================================================
% analyseResults.m
%==========================================================================
% Post-processing script for polymaxRIR batch results.
%
% PURPOSE
%   Compare two measurement conditions — an empty reverberant room and the
%   same room loaded with an absorption specimen — by analysing the modal
%   parameters extracted by polymaxRIR.m across all source–receiver
%   combinations in each condition.
%
% WORKFLOW
%   1.  Load all *_polymaxRIR.mat files from the two result folders.
%   2.  For each condition, pool poles from all IRs and vote on a COMMON
%       POLE BASIS: a frequency bin enters the basis only if it receives
%       a pole contribution from at least minFileFraction of the IRs.
%       The representative (f, c) for each basis pole is the median of the
%       contributing per-IR values.
%   3.  Re-fit least-squares modal amplitudes for every IR using its
%       condition's common basis (f and c are fixed; only a and b are
%       solved).  This places all IRs on a common parametrisation for
%       amplitude comparison.
%   4.  Compute per-display-band T60 statistics across all IRs for each
%       condition, from the original per-IR poleStats (not the re-fit).
%       Also extract the acoustic-indices T20 per IR for the same bands.
%   5.  Produce four figures:
%         Fig. 1 — Pole-frequency density histograms (LF and HF regimes)
%         Fig. 2 — Per-IR T60 spaghetti plots (inter-measurement spread)
%         Fig. 3 — MAIN RESULT: mean ± 1σ T60 per band, empty vs specimen,
%                  overlaid with acoustic-indices T20 for both conditions
%         Fig. 4 — Matched-basis analysis below fTransition:
%                  frequency shift and relative T60 change per mode
%   6.  Save all figures as .fig and one comparison .mat summary.
%
%
% DEPENDENCIES
%   Output .mat files from polymaxRIR.m  (batch mode, two folders).
%
% Author: Michele Ducceschi, March 2026
%==========================================================================

clear; close all; clc;

%% ========================== USER PARAMETERS ==============================

% --- Result folders -------------------------------------------------------
%   Set these to the subfolders created by polymaxRIR.m batch mode.
%   Each folder should contain one *_polymaxRIR.mat file per IR.
emptyDir    = './polymaxRIR_results/empty';
specimenDir = './polymaxRIR_results/St Gobain';

% --- Output ---------------------------------------------------------------
outDir  = './polymaxRIR_results/comparison';
verbose = true;

% --- Common-basis clustering ----------------------------------------------
%   binWidth_Hz     : frequency bin width for the voting histogram [Hz].
%                     Governs how close two per-IR pole estimates must be
%                     to be counted as the same physical mode.  Should be
%                     smaller than the expected mean mode spacing but larger
%                     than the pole-frequency estimation jitter across IRs
%                     (typically 0.5–2 Hz for reverberant rooms).
%   minFileFraction : fraction of IRs that must contribute a pole to a bin
%                     for that bin to enter the common basis.
%                     0.5 = majority vote (recommended starting point).
%                     Lower values yield a denser basis; higher values are
%                     more conservative and reject weakly excited modes.
binWidth_Hz     = 0.5;   % [Hz]
minFileFraction = 0.4;   % [0, 1]

% --- Display band resolution ----------------------------------------------
%   Controls the ISO band resolution used for ALL statistics, figures, and
%   the console table.  Independent of the LS fitting bands below.
%   'fullOctave'  — ISO 266 octave bands  (63, 125, 250 … Hz)
%   'thirdOctave' — ISO 266 third-octave bands
bandMode = 'thirdOctave';

% --- LS amplitude re-fit --------------------------------------------------
%   lsBandWidth_Hz  : width of the linear bands used for the LS amplitude
%                     solve [Hz].  Should be narrow enough that no single
%                     band contains too many basis poles — otherwise the
%                     normal-equation matrix becomes large and slow.
%                     Rule of thumb: similar to procBandWidth_Hz used in
%                     polymaxRIR (typically 25–100 Hz).  Wide octave display
%                     bands are NOT used for fitting — this parameter is
%                     independent of bandMode.
%   overlapFrac     : fractional extension of each LS window beyond the
%                     linear band edges (reduces boundary artefacts).
%   MinSecondsForFFT: FFT zero-padding floor [s].  Match the value used in
%                     polymaxRIR.m so that fv grids are identical.
lsBandWidth_Hz     = 50;    % [Hz]  linear LS band width
overlapFrac        = 0.10;
MinSecondsForFFT   = 10;

% --- Matched-pole comparison tolerance ------------------------------------
%   When pairing empty and specimen basis poles for Fig. 4, two poles are
%   matched if their frequency difference is within this tolerance [Hz].
matchTol_Hz = 3 * binWidth_Hz;

%% ========================== LOAD MAT FILES ===============================

% --- Check for existing synthesis results ---------------------------------
%   If comparison.mat already exists in outDir AND all per-IR
%   _commonBasis.mat files are present, skip the expensive pipeline and
%   load directly from disk.  Delete comparison.mat to force a full re-run.

compMatPath    = fullfile(outDir, 'comparison.mat');
emptySubDir    = fullfile(outDir, 'empty');
specimenSubDir = fullfile(outDir, 'specimen');

skipSynthesis = false;
if exist(compMatPath, 'file')
    fprintf('\nFound existing comparison.mat — checking per-IR cache...\n');
    nExpectedEmpty    = numel(dir(fullfile(emptyDir,    '*_polymaxRIR.mat')));
    nExpectedSpecimen = numel(dir(fullfile(specimenDir, '*_polymaxRIR.mat')));
    nCachedEmpty      = numel(dir(fullfile(emptySubDir,    '*_commonBasis.mat')));
    nCachedSpecimen   = numel(dir(fullfile(specimenSubDir, '*_commonBasis.mat')));
    if nCachedEmpty    == nExpectedEmpty    && nExpectedEmpty    > 0 && ...
       nCachedSpecimen == nExpectedSpecimen && nExpectedSpecimen > 0
        fprintf('Cache complete (%d + %d per-IR MATs) — skipping synthesis.\n', ...
            nCachedEmpty, nCachedSpecimen);
        skipSynthesis = true;
    else
        fprintf('Cache incomplete (%d/%d empty, %d/%d specimen) — running full pipeline.\n', ...
            nCachedEmpty, nExpectedEmpty, nCachedSpecimen, nExpectedSpecimen);
    end
end

if skipSynthesis
    fprintf('Loading comparison.mat...\n');
    compS = load(compMatPath);
    emptyBasis        = compS.emptyBasis;
    specimenBasis     = compS.specimenBasis;
    emptyFits         = compS.emptyFits;
    specimenFits      = compS.specimenFits;
    emptyT60mat       = compS.emptyT60mat;
    specimenT60mat    = compS.specimenT60mat;
    emptyMeanT60      = compS.emptyMeanT60;
    emptyStdT60       = compS.emptyStdT60;
    specimenMeanT60   = compS.specimenMeanT60;
    specimenStdT60    = compS.specimenStdT60;
    emptyMeanAIT20    = compS.emptyMeanAIT20;
    emptyStdAIT20     = compS.emptyStdAIT20;
    specimenMeanAIT20 = compS.specimenMeanAIT20;
    specimenStdAIT20  = compS.specimenStdAIT20;
    emptyBasisT60     = compS.emptyBasisT60;
    specimenBasisT60  = compS.specimenBasisT60;
    dispEdges         = compS.dispEdges;
    dispCentres       = compS.dispCentres;
    fTransition       = compS.fTransition;
    NDispBands        = numel(dispCentres);
    NEmpty            = numel(emptyFits.filename);
    NSpecimen         = numel(specimenFits.filename);
    fprintf('Loaded: %d empty | %d specimen | %d basis poles\n', ...
        NEmpty, NSpecimen, emptyBasis.NBasis);
else
    fprintf('\n=== LOADING RESULTS ===\n');
    emptyData    = loadConditionMats(emptyDir,    verbose);
    specimenData = loadConditionMats(specimenDir, verbose);

    NEmpty    = numel(emptyData);
    NSpecimen = numel(specimenData);

% Frequency range from the stored data
ref  = emptyData(1);
fLow = min(ref.dispEdgesVec);
fHigh = max(ref.dispEdgesVec);
fTransition = median([emptyData.fTransition]);

% --- Display bands (ISO octave or third-octave, rebuilt from bandMode) ----
[dispEdges, dispCentres] = buildISOBands(fLow, fHigh, bandMode);
NDispBands = numel(dispCentres);

% --- Linear LS bands (narrow fixed-width, used only for amplitude fitting) -
lsEdges  = buildLinearBands(fLow, fHigh, lsBandWidth_Hz);
NLSBands = numel(lsEdges) - 1;

fprintf('\nEmpty files    : %d\n', NEmpty);
fprintf('Specimen files : %d\n', NSpecimen);
fprintf('fTransition    : %.1f Hz  (median across empty files)\n', fTransition);
fprintf('Display bands  : %d  (%s)  |  [%.1f – %.1f] Hz\n', ...
    NDispBands, bandMode, dispEdges(1), dispEdges(end));
fprintf('LS bands       : %d linear bands of %.0f Hz  |  [%.1f – %.1f] Hz\n', ...
    NLSBands, lsBandWidth_Hz, lsEdges(1), lsEdges(end));

%% ==================== FIND COMMON POLE BASIS ============================
%
%   Modal frequencies are set by room geometry and are shared between both
%   conditions.  Only the decay times change when the specimen is added.
%   Strategy:
%     1. Vote on a common frequency basis from the EMPTY room IRs.
%     2. Re-estimate the damping for each basis frequency from the SPECIMEN
%        IRs (nearest pole within binWidth_Hz), keeping frequencies fixed.
%   Both conditions therefore share the same fBasis; only cBasis differs.

fprintf('\n=== FINDING COMMON POLE BASIS (from empty room) ===\n');
emptyBasis = findCommonBasis(emptyData, binWidth_Hz, minFileFraction, 'empty', verbose);

fprintf('\n=== ESTIMATING SPECIMEN DAMPING ON SHARED FREQUENCY BASIS ===\n');
specimenBasis = estimateDampingForBasis(emptyBasis, specimenData, binWidth_Hz, verbose);

%% ================ RE-FIT AMPLITUDES WITH COMMON BASIS ===================

fprintf('\n=== RE-FITTING AMPLITUDES (common basis, linear banded LS) ===\n');
emptyFits    = refitAmplitudes(emptyData,    emptyBasis,    lsEdges, overlapFrac, MinSecondsForFFT, verbose);
specimenFits = refitAmplitudes(specimenData, specimenBasis, lsEdges, overlapFrac, MinSecondsForFFT, verbose);

%% =================== PER-BAND T60 STATISTICS ============================

fprintf('\n=== PER-BAND T60 STATISTICS ===\n');

% NFiles × NDispBands matrix: entry (i,k) = mean T60 of all original poles
% in display band k for file i.  This captures the full per-IR distribution.
emptyT60mat    = extractT60matrix(emptyData,    dispEdges, NDispBands);
specimenT60mat = extractT60matrix(specimenData, dispEdges, NDispBands);

% Statistics across IRs (columns = bands)
emptyMeanT60    = mean(emptyT60mat,    1, 'omitnan');
emptyStdT60     = std( emptyT60mat,    0, 1, 'omitnan');
specimenMeanT60 = mean(specimenT60mat, 1, 'omitnan');
specimenStdT60  = std( specimenT60mat, 0, 1, 'omitnan');

% Acoustic-indices T20 per file, mapped to display band centres
emptyAIT20mat    = extractAIT20matrix(emptyData,    dispEdges, dispCentres, NDispBands);
specimenAIT20mat = extractAIT20matrix(specimenData, dispEdges, dispCentres, NDispBands);

emptyMeanAIT20    = mean(emptyAIT20mat,    1, 'omitnan');
emptyStdAIT20     = std( emptyAIT20mat,    0, 1, 'omitnan');
specimenMeanAIT20 = mean(specimenAIT20mat, 1, 'omitnan');
specimenStdAIT20  = std( specimenAIT20mat, 0, 1, 'omitnan');

% T60 derived directly from the common basis damping (single consensus value)
emptyBasisT60    = basisT60perBand(emptyBasis,    dispEdges, NDispBands);
specimenBasisT60 = basisT60perBand(specimenBasis, dispEdges, NDispBands);

if verbose
    fprintf('\n%8s  %9s  %9s  %9s  %9s  %9s  %9s\n', ...
        'fc(Hz)', 'E_mT60', 'E_sT60', 'S_mT60', 'S_sT60', 'E_AI', 'S_AI');
    for ib = 1 : NDispBands
        fprintf('%8.0f  %9.3f  %9.3f  %9.3f  %9.3f  %9.3f  %9.3f\n', ...
            dispCentres(ib), ...
            emptyMeanT60(ib),    emptyStdT60(ib), ...
            specimenMeanT60(ib), specimenStdT60(ib), ...
            emptyMeanAIT20(ib),  specimenMeanAIT20(ib));
    end
end

end   % if ~skipSynthesis

%% ============================= FIGURES ==================================

colEmpty    = [0.18 0.40 0.75];   % blue
colSpecimen = [0.82 0.20 0.15];   % red
colSch      = [0.10 0.55 0.10];   % green

%--------------------------------------------------------------------------
%  Figure 1 — T60 distribution comparison (LF and HF)
%
%   Frequency histograms are uninformative because both conditions share
%   the same frequency basis.  The meaningful comparison is the T60
%   distribution: frequencies are geometry-locked, damping is absorption-driven.
%--------------------------------------------------------------------------
hFig1 = figure('Name', 'T60 distributions', ...
    'Color', 'w', 'Position', [50 50 1300 520]);

% T60 for every basis pole in each condition
T60_E_all = 3*log(10) ./ max(emptyBasis.cBasis,    1e-10);
T60_S_all = 3*log(10) ./ max(specimenBasis.cBasis,  1e-10);

% Split by regime
lfMask = emptyBasis.fBasis < fTransition;   % same for both (shared fBasis)
hfMask = ~lfMask;

T60_E_LF = T60_E_all(lfMask);   T60_S_LF = T60_S_all(lfMask);
T60_E_HF = T60_E_all(hfMask);   T60_S_HF = T60_S_all(hfMask);

% Shared bin edges per regime (same for both conditions so they overlay cleanly)
T60binLF = linspace(0, max([T60_E_LF; T60_S_LF]) * 1.05, 40);
T60binHF = linspace(0, max([T60_E_HF; T60_S_HF]) * 1.05, 40);

subplot(1, 2, 1);
hold on;
histogram(T60_E_LF, T60binLF, 'FaceColor', colEmpty,    'EdgeColor', 'none', ...
    'Normalization', 'probability', 'DisplayName', 'Empty');
histogram(T60_S_LF, T60binLF, 'FaceColor', colSpecimen, 'EdgeColor', 'none', ...
    'FaceAlpha', 0.6, 'Normalization', 'probability', 'DisplayName', 'Specimen');
xline(median(T60_E_LF), '--', 'Color', colEmpty,    'LineWidth', 1.5, ...
    'Label', sprintf('E med=%.2fs', median(T60_E_LF)), 'HandleVisibility', 'off');
xline(median(T60_S_LF), '--', 'Color', colSpecimen, 'LineWidth', 1.5, ...
    'Label', sprintf('S med=%.2fs', median(T60_S_LF)), 'HandleVisibility', 'off');
xlabel('T_{60} (s)');  ylabel('Probability');
legend('show', 'Location', 'northeast');
title(sprintf('LF regime  (f < %.0f Hz)  —  %d shared modes', fTransition, nnz(lfMask)));
grid on;  box on;

subplot(1, 2, 2);
hold on;
histogram(T60_E_HF, T60binHF, 'FaceColor', colEmpty,    'EdgeColor', 'none', ...
    'Normalization', 'probability', 'DisplayName', 'Empty');
histogram(T60_S_HF, T60binHF, 'FaceColor', colSpecimen, 'EdgeColor', 'none', ...
    'FaceAlpha', 0.6, 'Normalization', 'probability', 'DisplayName', 'Specimen');
xline(median(T60_E_HF), '--', 'Color', colEmpty,    'LineWidth', 1.5, ...
    'Label', sprintf('E med=%.2fs', median(T60_E_HF)), 'HandleVisibility', 'off');
xline(median(T60_S_HF), '--', 'Color', colSpecimen, 'LineWidth', 1.5, ...
    'Label', sprintf('S med=%.2fs', median(T60_S_HF)), 'HandleVisibility', 'off');
xlabel('T_{60} (s)');  ylabel('Probability');
legend('show', 'Location', 'northeast');
title(sprintf('HF regime  (f \\geq %.0f Hz)  —  %d shared modes', fTransition, nnz(hfMask)));
grid on;  box on;

sgtitle('T_{60} distribution: empty vs specimen  (shared frequency basis)', 'FontWeight', 'bold');

%--------------------------------------------------------------------------
%  Figure 2 — Per-IR T60 spaghetti: both conditions in one panel
%--------------------------------------------------------------------------
hFig2 = figure('Name', 'T60 spaghetti', ...
    'Color', 'w', 'Position', [50 620 900 500]);

hold on;
for iF = 1 : NEmpty
    h = semilogx(dispCentres, emptyT60mat(iF, :), '-', ...
        'Color', colEmpty, 'LineWidth', 0.8, 'HandleVisibility', 'off');
    h.Color(4) = 0.20;
end
for iF = 1 : NSpecimen
    h = semilogx(dispCentres, specimenT60mat(iF, :), '-', ...
        'Color', colSpecimen, 'LineWidth', 0.8, 'HandleVisibility', 'off');
    h.Color(4) = 0.20;
end
semilogx(dispCentres, emptyMeanT60, '-', ...
    'Color', colEmpty,    'LineWidth', 2.5, 'DisplayName', 'Empty — mean');
semilogx(dispCentres, specimenMeanT60, '-', ...
    'Color', colSpecimen, 'LineWidth', 2.5, 'DisplayName', 'Specimen — mean');
xline(fTransition, '--', 'Color', colSch, 'LineWidth', 1.2, 'Label', 'f_{Sch}');
xlabel('f_c (Hz)');  ylabel('T_{60} (s)');
legend('show', 'Location', 'best');
set(gca, 'XScale', 'log');  grid on;  box on;
xlim([dispEdges(1), dispEdges(end)]);
title('Inter-IR T_{60} variability (faint) vs mean (solid) — empty and specimen', ...
    'FontWeight', 'bold');

%--------------------------------------------------------------------------
%  Figure 3 — MAIN RESULT: mean ± 1σ T60, both conditions + AI T20
%--------------------------------------------------------------------------
hFig3 = figure('Name', 'T60 comparison — main result', ...
    'Color', 'w', 'Position', [1400 50 900 620]);

hold on;

% Shade the LF region (below fTransition)
yLim_guess = max([emptyMeanT60, specimenMeanT60, ...
    emptyMeanAIT20, specimenMeanAIT20], [], 'omitnan') * 1.25;
if isnan(yLim_guess) || yLim_guess <= 0, yLim_guess = 5; end
patch([dispEdges(1), fTransition, fTransition, dispEdges(1)], ...
    [0, 0, yLim_guess*1.1, yLim_guess*1.1], ...
    [0.92 0.96 0.92], 'EdgeColor', 'none', 'FaceAlpha', 0.6, ...
    'HandleVisibility', 'off');

% --- Empty room: modal T60 mean ± 1σ ---
eb_E = errorbar(dispCentres, emptyMeanT60, emptyStdT60, ...
    'o-', ...
    'Color',           colEmpty, ...
    'MarkerFaceColor', colEmpty, ...
    'MarkerSize',      7, ...
    'LineWidth',       1.8, ...
    'CapSize',         7, ...
    'DisplayName',     'Empty — modal T_{60}  (mean \pm 1\sigma)');

% --- Specimen room: modal T60 mean ± 1σ ---
eb_S = errorbar(dispCentres, specimenMeanT60, specimenStdT60, ...
    's-', ...
    'Color',           colSpecimen, ...
    'MarkerFaceColor', colSpecimen, ...
    'MarkerSize',      7, ...
    'LineWidth',       1.8, ...
    'CapSize',         7, ...
    'DisplayName',     'Specimen — modal T_{60}  (mean \pm 1\sigma)');

% --- Empty room: acoustic-indices T20 mean ± 1σ ---
errorbar(dispCentres, emptyMeanAIT20, emptyStdAIT20, ...
    'o--', ...
    'Color',           colEmpty * 0.65, ...
    'MarkerFaceColor', 'none', ...
    'MarkerSize',      7, ...
    'LineWidth',       1.2, ...
    'CapSize',         5, ...
    'DisplayName',     'Empty — AI T_{20}  (mean \pm 1\sigma)');

% --- Specimen room: acoustic-indices T20 mean ± 1σ ---
errorbar(dispCentres, specimenMeanAIT20, specimenStdAIT20, ...
    's--', ...
    'Color',           colSpecimen * 0.65, ...
    'MarkerFaceColor', 'none', ...
    'MarkerSize',      7, ...
    'LineWidth',       1.2, ...
    'CapSize',         5, ...
    'DisplayName',     'Specimen — AI T_{20}  (mean \pm 1\sigma)');

% --- Schroeder transition ---
xline(fTransition, '-', 'Color', colSch, 'LineWidth', 1.5, ...
    'Label', 'f_{Sch}', 'LabelVerticalAlignment', 'top');

set(gca, 'XScale', 'log');
xlabel('Band centre frequency  f_c (Hz)', 'FontSize', 12);
ylabel('T_{60}  (s)',                     'FontSize', 12);
xlim([dispEdges(1), dispEdges(end)]);
ylim([0, yLim_guess * 1.1]);
legend('show', 'Location', 'best', 'FontSize', 10);
grid on;  box on;

text(fTransition * 0.55, yLim_guess * 0.96, 'LF  (modal)', ...
    'HorizontalAlignment', 'right', 'Color', [0.3 0.55 0.3], 'FontSize', 9);
text(fTransition * 1.45, yLim_guess * 0.96, 'HF  (diffuse)', ...
    'HorizontalAlignment', 'left',  'Color', [0.3 0.55 0.3], 'FontSize', 9);

title({'T_{60} per band: empty vs specimen', ...
    '(solid = PolyMax modal fit,  dashed = acoustic-indices T_{20})'}, ...
    'FontSize', 12);

%% ============================== SAVE =====================================

if ~exist(outDir,        'dir'), mkdir(outDir);        end
if ~exist(emptySubDir,   'dir'), mkdir(emptySubDir);   end
if ~exist(specimenSubDir,'dir'), mkdir(specimenSubDir);end

% --- Figures (always regenerated) -----------------------------------------
figNames   = {'fig1_T60distributions.fig', 'fig2_T60spaghetti.fig', ...
              'fig3_T60comparison.fig'};
figHandles = {hFig1, hFig2, hFig3};
for iF = 1 : 3
    savefig(figHandles{iF}, fullfile(outDir, figNames{iF}));
end
fprintf('\nFigures saved:\n');
for iF = 1 : 3
    fprintf('  %s\n', fullfile(outDir, figNames{iF}));
end

% --- Per-IR MAT and WAV (only on full pipeline run) -----------------------
if ~skipSynthesis
    fprintf('\nPer-IR outputs:\n');
    for iF = 1 : numel(emptyFits.filename)
        fname = emptyFits.filename{iF};  fs_ir = emptyFits.fs(iF);
        irS.filename = fname;  irS.fBasis = emptyBasis.fBasis;
        irS.cBasis   = emptyBasis.cBasis;
        irS.aBasis   = emptyFits.aBasis{iF};  irS.bBasis = emptyFits.bBasis{iF};
        irS.nRMSE    = emptyFits.nRMSE(iF);   irS.fs     = fs_ir;
        irS.condition = 'empty';
        save(fullfile(emptySubDir,    [fname, '_commonBasis.mat']), '-struct', 'irS');
        audiowrite(fullfile(emptySubDir, [fname, '_commonBasis_SYNTH.wav']), emptyFits.outSyn{iF}, fs_ir);
        fprintf('  [empty]    %s\n', fname);
    end
    for iF = 1 : numel(specimenFits.filename)
        fname = specimenFits.filename{iF};  fs_ir = specimenFits.fs(iF);
        irS.filename = fname;  irS.fBasis = specimenBasis.fBasis;
        irS.cBasis   = specimenBasis.cBasis;
        irS.aBasis   = specimenFits.aBasis{iF};  irS.bBasis = specimenFits.bBasis{iF};
        irS.nRMSE    = specimenFits.nRMSE(iF);   irS.fs     = fs_ir;
        irS.condition = 'specimen';
        save(fullfile(specimenSubDir, [fname, '_commonBasis.mat']), '-struct', 'irS');
        audiowrite(fullfile(specimenSubDir, [fname, '_commonBasis_SYNTH.wav']), specimenFits.outSyn{iF}, fs_ir);
        fprintf('  [specimen] %s\n', fname);
    end

    % --- comparison.mat (only on full pipeline run) -----------------------
    compS.emptyBasis        = emptyBasis;     compS.specimenBasis     = specimenBasis;
    compS.emptyFits         = emptyFits;      compS.specimenFits      = specimenFits;
    compS.emptyT60mat       = emptyT60mat;    compS.specimenT60mat    = specimenT60mat;
    compS.emptyMeanT60      = emptyMeanT60;   compS.emptyStdT60       = emptyStdT60;
    compS.specimenMeanT60   = specimenMeanT60;compS.specimenStdT60    = specimenStdT60;
    compS.emptyMeanAIT20    = emptyMeanAIT20; compS.emptyStdAIT20     = emptyStdAIT20;
    compS.specimenMeanAIT20 = specimenMeanAIT20; compS.specimenStdAIT20 = specimenStdAIT20;
    compS.emptyBasisT60     = emptyBasisT60;  compS.specimenBasisT60  = specimenBasisT60;
    compS.dispEdges         = dispEdges;      compS.dispCentres       = dispCentres;
    compS.fTransition       = fTransition;
    compS.emptyDir          = emptyDir;       compS.specimenDir       = specimenDir;
    compS.bandMode          = bandMode;       compS.lsBandWidth_Hz    = lsBandWidth_Hz;
    compS.binWidth_Hz       = binWidth_Hz;    compS.minFileFraction   = minFileFraction;
    compS.processedDate     = datestr(now);   %#ok<TNOW1,DATST>
    save(compMatPath, '-struct', 'compS');
    fprintf('\ncomparison.mat saved : %s\n', compMatPath);
else
    fprintf('\nPer-IR files and comparison.mat not overwritten (loaded from cache).\n');
end

fprintf('\n=== DONE ===\n');
if skipSynthesis
    fprintf('Fast path: synthesis skipped, figures regenerated.\n');
else
    fprintf('Full pipeline: %d empty + %d specimen IRs processed.\n', ...
        numel(emptyFits.filename), numel(specimenFits.filename));
end
fprintf('Output folder : %s\n', outDir);


%% ########################################################################
%%  LOCAL FUNCTIONS
%% ########################################################################

% =========================================================================
%  Data loading
% =========================================================================

function data = loadConditionMats(folder, verbose)
%LOADCONDITIONMATS  Load all *_polymaxRIR.mat files from FOLDER.

matFiles = dir(fullfile(folder, '*_polymaxRIR.mat'));
if isempty(matFiles)
    error('analyseResults:noMats', ...
        'No *_polymaxRIR.mat files found in ''%s''.', folder);
end
N = numel(matFiles);
if verbose
    fprintf('  Loading %d file(s) from %s\n', N, folder);
end
for k = 1 : N
    data(k) = load(fullfile(folder, matFiles(k).name)); %#ok<AGROW>
    if verbose
        fprintf('    %2d  %s\n', k, matFiles(k).name);
    end
end
end

% =========================================================================
%  Common pole basis
% =========================================================================

function basis = findCommonBasis(dataArray, binWidth, minFrac, label, verbose)
%FINDCOMMONBASIS  Vote-based frequency clustering across all IRs.
%
%   For each frequency bin of width BINWIDTH, counts how many IRs contribute
%   at least one pole.  Bins with contributions from >= MINFRAC of IRs enter
%   the basis; their representative (f, c) is the median of the per-IR values
%   closest to the bin centre.

NFiles   = numel(dataArray);
minCount = ceil(minFrac * NFiles);

fLo = min(dataArray(1).dispEdgesVec);
fHi = max(dataArray(1).dispEdgesVec);
binEdges = fLo : binWidth : fHi + binWidth;
NBins    = numel(binEdges) - 1;

fBuf = zeros(NBins, 1);
cBuf = zeros(NBins, 1);
valid = false(NBins, 1);

for ib = 1 : NBins
    bLo  = binEdges(ib);
    bHi  = binEdges(ib + 1);
    bCtr = 0.5 * (bLo + bHi);

    fRep  = nan(NFiles, 1);
    cRep  = nan(NFiles, 1);
    count = 0;

    for k = 1 : NFiles
        inBin = (dataArray(k).fAll >= bLo) & (dataArray(k).fAll < bHi);
        if ~any(inBin), continue; end
        count = count + 1;
        fThis = dataArray(k).fAll(inBin);
        cThis = dataArray(k).cAll(inBin);
        [~, iBest] = min(abs(fThis - bCtr));
        fRep(count) = fThis(iBest);
        cRep(count) = cThis(iBest);
    end

    if count >= minCount
        valid(ib)  = true;
        fBuf(ib)   = median(fRep(1:count), 'omitnan');
        cBuf(ib)   = median(cRep(1:count), 'omitnan');
    end
end

basis.fBasis   = fBuf(valid);
basis.cBasis   = cBuf(valid);
basis.NBasis   = nnz(valid);
basis.label    = label;
basis.binWidth = binWidth;
basis.minFrac  = minFrac;

if verbose
    fprintf('  [%s]  basis poles = %d  (threshold: >= %d/%d files per bin)\n', ...
        label, basis.NBasis, minCount, NFiles);
end
end

% =========================================================================

function specimenBasis = estimateDampingForBasis(emptyBasis, specimenData, binWidth, verbose)
%ESTIMATEDAMPINGFORBASIS  Re-estimate damping on a fixed frequency basis.
%
%   Takes the shared frequency basis established from the empty room and
%   estimates per-basis-pole damping coefficients from the specimen IRs.
%   For each basis frequency f_k, searches each specimen IR for the nearest
%   pole within BINWIDTH Hz and takes the median c across all IRs that
%   contribute.  Basis poles with no specimen coverage retain the empty-room
%   damping as a fallback (flagged in the console).
%
%   The returned specimenBasis has the same fBasis as emptyBasis — frequencies
%   are NOT re-voted — only cBasis differs.

NFiles  = numel(specimenData);
NBasis  = emptyBasis.NBasis;
fBasis  = emptyBasis.fBasis;   % shared, fixed

cMat   = nan(NBasis, NFiles);   % per-pole per-file damping estimates

for iF = 1 : NFiles
    fAll = specimenData(iF).fAll(:);
    cAll = specimenData(iF).cAll(:);
    for ip = 1 : NBasis
        [dMin, idx] = min(abs(fAll - fBasis(ip)));
        if dMin <= binWidth
            cMat(ip, iF) = cAll(idx);
        end
    end
end

% Median across IRs; fall back to empty-room damping where coverage is missing
cBasisNew  = median(cMat, 2, 'omitnan');
nFallback  = nnz(isnan(cBasisNew));
if nFallback > 0
    if verbose
        fprintf('  [specimen]  %d/%d basis poles have no specimen coverage — using empty-room damping\n', ...
            nFallback, NBasis);
    end
    cBasisNew(isnan(cBasisNew)) = emptyBasis.cBasis(isnan(cBasisNew));
end

specimenBasis          = emptyBasis;   % copy structure (same fBasis, metadata)
specimenBasis.cBasis   = cBasisNew;
specimenBasis.label    = 'specimen';

if verbose
    T60empty    = 3*log(10) ./ max(emptyBasis.cBasis,   1e-10);
    T60specimen = 3*log(10) ./ max(specimenBasis.cBasis, 1e-10);
    fprintf('  [specimen]  median T60 change: %.3f s -> %.3f s  (%.0f%% reduction)\n', ...
        median(T60empty, 'omitnan'), median(T60specimen, 'omitnan'), ...
        100 * (1 - median(T60specimen,'omitnan') / median(T60empty,'omitnan')));
end
end

% =========================================================================
%  Amplitude re-fit with common basis
% =========================================================================

function fits = refitAmplitudes(dataArray, basis, lsEdges, overlapFrac, minSecsFFT, verbose)
%REFITAMPLITUDES  Banded least-squares amplitude solve with a fixed (f,c) basis.
%
%   For each IR, recomputes the FFT from the stored hAligned, then solves for
%   residue amplitudes (a, b) of the common-basis poles using narrow LINEAR
%   bands of width lsBandWidth_Hz (with OVERLAPFRAC overlap on each side).
%   Linear bands are used — NOT octave display bands — so that no single LS
%   solve involves an unmanageably large number of poles.

NFiles  = numel(dataArray);
NBands  = numel(lsEdges) - 1;
fBasis  = basis.fBasis;
cBasis  = basis.cBasis;
NBasis  = basis.NBasis;

% Initialise output as scalar struct with array fields (avoids struct-array trap)
fits.filename = cell(NFiles, 1);
fits.aBasis   = cell(NFiles, 1);
fits.bBasis   = cell(NFiles, 1);
fits.outSyn   = cell(NFiles, 1);
fits.fs       = nan(NFiles, 1);
fits.nRMSE    = nan(NFiles, 1);

for iF = 1 : NFiles
    fname = dataArray(iF).filename;
    fprintf('  [%s]  IR %d/%d : %s\n', basis.label, iF, NFiles, fname);

    h    = dataArray(iF).refAI.hAligned(:);
    fs   = dataArray(iF).fs;
    Lfft = max(round(minSecsFFT * fs), numel(h));
    [spc, fv, ~] = takefft(h, fs, Lfft);

    aBasis  = zeros(NBasis, 1);
    bBasis  = zeros(NBasis, 1);
    spcSyn  = zeros(size(spc));   % accumulate band contributions here

    for ib = 1 : NBands
        flow  = lsEdges(ib);
        fhigh = lsEdges(ib + 1);
        bw    = fhigh - flow;

        inBand = (fBasis >= flow) & (fBasis < fhigh);
        Np = nnz(inBand);
        if Np == 0, continue; end

        fLo_w = max(fv(2),   flow  - overlapFrac * bw);
        fHi_w = min(fv(end), fhigh + overlapFrac * bw);
        [~, is] = min(abs(fv - fLo_w));
        [~, ie] = min(abs(fv - fHi_w));
        ie = max(ie, is + 1);

        omvSlice = 2*pi * fv(is:ie);
        spcSlice = spc(is:ie);

        [aB, bB] = leastSqWeights(omvSlice, Np, ...
            fBasis(inBand), cBasis(inBand), spcSlice);

        idx = find(inBand);
        aBasis(idx) = aB;
        bBasis(idx) = bB;

        % Synthesise this band's contribution into the full spectrum vector.
        % Only the poles in this band are used, over the core (non-overlapped)
        % bin range so that adjacent bands don't double-count.
        [~, is_core] = min(abs(fv - flow));
        [~, ie_core] = min(abs(fv - fhigh));
        ie_core = max(ie_core, is_core + 1);
        omvCore = 2*pi * fv(is_core:ie_core);
        spcSyn(is_core:ie_core) = spcSyn(is_core:ie_core) + ...
            spectrumBuildBand(aB, bB, cBasis(inBand), fBasis(inBand), omvCore);
    end

    fits.filename{iF} = dataArray(iF).filename;
    fits.aBasis{iF}   = aBasis;
    fits.bBasis{iF}   = bBasis;
    fits.fs(iF)       = fs;

    % nRMSE computed from the band-assembled synthetic spectrum —
    % no full spectrumBuild call needed.
    outRef  = real(ifft([spc; 0; conj(spc(end:-1:2))])) * fs;
    Lout    = numel(outRef);
    outSyn  = synthesizeIR(aBasis, bBasis, basis.cBasis, basis.fBasis, fs, Lout);
    fits.nRMSE(iF)   = normalisedRMSE_time(outRef, outSyn);
    fits.outSyn{iF}  = outSyn / (max(abs(outSyn)) + eps);
    fprintf('    -> nRMSE (common basis) = %.4g\n', fits.nRMSE(iF));
end
end

% =========================================================================
%  Per-band T60 extraction
% =========================================================================

function T60mat = extractT60matrix(dataArray, dispEdges, NDispBands)
%EXTRACTT60MATRIX  NFiles × NDispBands matrix of per-band mean T60.
%
%   Entry (i,k) is the mean T60 of all poles from file i whose frequency
%   falls in display band k.  NaN when a band has no poles in that file.

NFiles  = numel(dataArray);
T60mat  = nan(NFiles, NDispBands);

for iF = 1 : NFiles
    fAllF = dataArray(iF).fAll(:);
    cAllF = dataArray(iF).cAll(:);
    for ib = 1 : NDispBands
        inBand = (fAllF >= dispEdges(ib)) & (fAllF < dispEdges(ib + 1));
        if ~any(inBand), continue; end
        cBand = cAllF(inBand);
        cBand = cBand(cBand > 0);
        if isempty(cBand), continue; end
        T60mat(iF, ib) = mean(3*log(10) ./ cBand, 'omitnan');
    end
end
end

% -------------------------------------------------------------------------

function AIT20mat = extractAIT20matrix(dataArray, dispEdges, dispCentres, NDispBands)
%EXTRACTAIT20MATRIX  NFiles × NDispBands matrix of acoustic-indices T20.
%
%   Maps each file's refAI.T20 values to the display bands by matching
%   AI band centres to display band centres in log-frequency.

NFiles   = numel(dataArray);
AIT20mat = nan(NFiles, NDispBands);

for iF = 1 : NFiles
    ai_fc  = dataArray(iF).refAI.bandCentersHz(:);
    ai_T20 = dataArray(iF).refAI.T20(:);
    for ib = 1 : NDispBands
        [~, idx] = min(abs(log2(ai_fc / max(dispCentres(ib), eps))));
        if ~isempty(idx) && ~isnan(ai_T20(idx)) && ai_T20(idx) > 0
            AIT20mat(iF, ib) = ai_T20(idx);
        end
    end
end
end

% -------------------------------------------------------------------------

function basisT60 = basisT60perBand(basis, dispEdges, NDispBands)
%BASIST60PERBAND  Mean T60 of common-basis poles per display band.

basisT60 = nan(1, NDispBands);
for ib = 1 : NDispBands
    inBand = (basis.fBasis >= dispEdges(ib)) & ...
             (basis.fBasis <  dispEdges(ib + 1));
    if ~any(inBand), continue; end
    c = basis.cBasis(inBand);
    c = c(c > 0);
    if isempty(c), continue; end
    basisT60(ib) = mean(3*log(10) ./ c, 'omitnan');
end
end

% =========================================================================
%  Matched-basis comparison
% =========================================================================

function [fE_match, cE_match, fS_match, cS_match] = ...
    matchBasisPoles(emptyBasis, specimenBasis, fTransition, tolHz)
%MATCHBASISPOLES  Pair LF basis poles between conditions by proximity.
%
%   For each empty-room basis pole below FTRANSITION, finds the nearest
%   specimen-room basis pole within TOLHZ.  Unmatched poles are discarded.

fE = emptyBasis.fBasis;
cE = emptyBasis.cBasis;
fS = specimenBasis.fBasis;
cS = specimenBasis.cBasis;

% Restrict to LF regime
lfE = fE < fTransition;
fE  = fE(lfE);  cE = cE(lfE);

fE_match = [];  cE_match = [];
fS_match = [];  cS_match = [];

for ip = 1 : numel(fE)
    [dMin, idx] = min(abs(fS - fE(ip)));
    if dMin <= tolHz
        fE_match(end+1) = fE(ip);           %#ok<AGROW>
        cE_match(end+1) = cE(ip);           %#ok<AGROW>
        fS_match(end+1) = fS(idx);          %#ok<AGROW>
        cS_match(end+1) = cS(idx);          %#ok<AGROW>
    end
end

fE_match = fE_match(:);  cE_match = cE_match(:);
fS_match = fS_match(:);  cS_match = cS_match(:);
end

% =========================================================================
%  Signal processing utilities  (copied from polymaxRIR.m)
% =========================================================================

function [edges, fc] = buildISOBands(fLow, fHigh, bandMode)
%BUILDISOBANDS  ISO 266 octave or third-octave band edges for FRANGE.
%
%   Bands are selected by centre frequency falling within [fLow, fHigh].
%   Edges are NOT clipped so that every band has standard ISO width.

switch lower(bandMode)
    case 'fulloctave',   B = 1;
    case 'thirdoctave',  B = 3;
    otherwise
        error('analyseResults:unknownBandMode', ...
            'bandMode must be ''fullOctave'' or ''thirdOctave''.');
end

halfStep = 2^(1 / (2*B));
nMin = floor(B * log2(fLow  / 1000)) - 1;
nMax = ceil( B * log2(fHigh / 1000)) + 1;
n    = nMin : nMax;
fc_all = 1000 * 2.^(n / B);
mask   = (fc_all >= fLow) & (fc_all <= fHigh);
fc     = fc_all(mask);
lo     = fc / halfStep;
hi     = fc * halfStep;
edges  = [lo(:).', hi(end)];
end

% -------------------------------------------------------------------------

function edges = buildLinearBands(fLow, fHigh, bwHz)
%BUILDLINEARBANDS  Narrow linear bands of fixed width BWHZ over [fLow, fHigh].

edges = fLow : bwHz : fHigh;
if edges(end) < fHigh
    edges(end+1) = fHigh;
end
edges = edges(:).';
end

% -------------------------------------------------------------------------

function h = synthesizeIR(aVec, bVec, cVec, fVec, fs, Nout)
%SYNTHESIZEIR  Direct time-domain modal synthesis — no IFFT wrap-around.
%
%   h(n) = 2·Re[ Σ_p (a_p+j·b_p) · exp(λ_p · n/fs) ]
%   λ_p  = -c_p + j·√((2πf_p)² - c_p²)
%
%   IFFT-based reconstruction aliases energy from long-decay poles (T60
%   exceeding the FFT window) back onto the end of the block.  This function
%   produces the correct causal decay to arbitrary length.

aVec = aVec(:);  bVec = bVec(:);  cVec = cVec(:);  fVec = fVec(:);
Np     = numel(aVec);
rad    = max((2*pi*fVec).^2 - cVec.^2, 1e-12);
lambda = -cVec + 1j*sqrt(rad);
resids = aVec + 1j*bVec;
n      = (0 : Nout-1);
chunkSize = max(1, floor(50e6 / (Nout * 16)));
h = zeros(Nout, 1);
for p = 1 : chunkSize : Np
    idx = p : min(p + chunkSize - 1, Np);
    h   = h + 2 * real(sum(resids(idx) .* exp(lambda(idx) .* n / fs), 1)).';
end
end

% -------------------------------------------------------------------------

function [aBand, bBand] = leastSqWeights(omvBand, Npoles, fBand, cBand, fftRef)
%LEASTSQWEIGHTS  Least-squares residue amplitudes for a modal partial-fraction model.
Lfv   = numel(omvBand);
fBand = fBand(:).';
cBand = cBand(:).';
rad   = max((2*pi*fBand).^2 - cBand.^2, 1e-12);
poles = -cBand + 1j*sqrt(rad);
jomv  = 1j * omvBand(:);
D     = jomv - poles;
Dc    = jomv - conj(poles);
Ra    = 1./D  + 1./Dc;
Rb    = 1j./D - 1j./Dc;
R     = [Ra, Rb];
RtR   = R' * R;
eps_rel = 1e-8;
reg   = eps_rel * trace(RtR) / (2*Npoles) * eye(2*Npoles);
solb  = (RtR + reg) \ (R' * fftRef(:));
aT    = solb(1:Npoles);
bT    = solb(Npoles + 1 : end);
aBand = real(aT) - imag(bT);
bBand = real(bT) + imag(aT);
end

% -------------------------------------------------------------------------

function spc = spectrumBuild(aVec, bVec, cVec, fVec, omv)
%SPECTRUMBUILD  Chunked modal partial-fraction spectrum evaluation.
aVec   = aVec(:);  bVec = bVec(:);  cVec = cVec(:);  fVec = fVec(:);
omv    = omv(:).';
Np     = numel(aVec);
Nfft   = numel(omv);
rad    = max((2*pi*fVec).^2 - cVec.^2, 1e-12);
poles  = -cVec + 1j*sqrt(rad);
resids = aVec + 1j*bVec;
chunkSize = max(1, floor(50e6 / (Nfft * 16)));
spc = zeros(1, Nfft);
for p = 1 : chunkSize : Np
    idx = p : min(p + chunkSize - 1, Np);
    r   = resids(idx);
    pl  = poles(idx);
    D   = 1j*omv - pl;
    Dc  = 1j*omv - conj(pl);
    spc = spc + sum(r./D + conj(r)./Dc, 1);
end
spc = spc(:);
end

% -------------------------------------------------------------------------

function nrmse = normalisedRMSE_time(xRef, xSyn)
%NORMALISEDRMSE_TIME  nRMSE = RMS(xSyn-xRef) / max(|xRef|).
xRef  = xRef(:);  xSyn = xSyn(:);
L     = min(numel(xRef), numel(xSyn));
nrmse = sqrt(mean((xSyn(1:L) - xRef(1:L)).^2)) / (max(abs(xRef(1:L))) + eps);
end

% -------------------------------------------------------------------------

function spc = spectrumBuildBand(aVec, bVec, cVec, fVec, omv)
%SPECTRUMBUILDBBAND  Modal spectrum for a single band (small Np, no chunking needed).
%
%   Identical algebra to spectrumBuild but operates on a short omv slice
%   with only the Np poles belonging to one display band, so the Np×Nbins
%   broadcast matrix is always small and memory-safe.

aVec = aVec(:);  bVec = bVec(:);  cVec = cVec(:);  fVec = fVec(:);
omv  = omv(:).';
rad    = max((2*pi*fVec).^2 - cVec.^2, 1e-12);
poles  = -cVec + 1j*sqrt(rad);   % Np×1
resids = aVec  + 1j*bVec;        % Np×1
D   = 1j*omv - poles;            % Np×Nbins
Dc  = 1j*omv - conj(poles);
spc = sum(resids./D + conj(resids)./Dc, 1);
spc = spc(:);
end

% -------------------------------------------------------------------------

function [fftinput, fv, omv] = takefft(input, fs, fftlength)
%TAKEFFT  Single-sided FFT normalised by 1/fs.
N        = max(fftlength, numel(input));
fftinput = fft(input, N) / fs;
fftinput = fftinput(1 : floor(N/2));
fv       = (0 : N-1).' * fs / N;
fv       = fv(1 : floor(N/2));
omv      = 2*pi*fv;
end