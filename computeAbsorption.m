%==========================================================================
% computeAbsorption.m
%==========================================================================
% Sound absorption coefficient of a specimen from polymaxRIR comparison
% results, combining per-mode modal estimates (LF) and Sabine band estimates
% (HF) into a single wideband α curve.
%
% THEORY
%   In Sabine's room-acoustic model every mode p decays with rate:
%
%       c_p  =  c_sound · A_total / (8 V)
%
%   where A_total = Σ_i α_i S_i is the total absorption area [m²].
%   Adding a specimen of area S_spec and coefficient α_spec raises A_total
%   by ΔA = α_spec · S_spec, so:
%
%       Δc_p  =  c_S,p  −  c_E,p  =  c_sound · α_spec · S_spec / (8 V)
%
%   Rearranging gives the per-mode estimate:
%
%       α_spec,p  =  8 V · Δc_p / (c_sound · S_spec)            [LF method]
%
%   For HF bands (above fTransition) individual modes are unresolved.
%   The standard Sabine formula applied to band-averaged T60 gives:
%
%       A_total  =  55.3 V / (c_sound · T60)
%       ΔA       =  55.3 V / c_sound · (1/T60_S − 1/T60_E)
%       α_spec   =  ΔA / S_spec                                  [HF method]
%
% ADVANTAGES OVER ISO SCHROEDER / SABINE METHOD
%   1. Direct use of c (linear damping) avoids the nonlinearity bias from
%      averaging T60 and then inverting; c is the physically fundamental
%      quantity that Sabine's model is written in.
%   2. Per-mode weighting by residue amplitude |r_p| = √(a_p²+b_p²):
%      well-excited modes (large |r|) contribute more; poorly resolved ones
%      are down-weighted automatically.
%   3. Uncertainty propagated from the inter-IR spread in cBasis:
%      the standard deviation of α across modes in a band gives a direct
%      confidence measure, not a post-hoc curve-fit error.
%   4. Specimen absorption extracted from Δc — wall absorption cancels
%      exactly without knowing α_walls; only ΔA = α_spec · S_spec remains.
%   5. No noise floor dependency: Schroeder integration fails when the
%      late-decay SNR is low; Δc is estimated from the stability diagram,
%      which is robust to background noise in the tails.
%
% INPUTS (from comparison.mat produced by analyseResults.m)
%   Plus user-supplied room geometry (V, c_sound, S_spec).
%
% OUTPUTS
%   α per display band (LF: modal, HF: Sabine), with uncertainty bounds.
%   Saved to outDir/absorption.mat and outDir/absorption_alpha.fig.
%
% DEPENDENCIES
%   comparison.mat from analyseResults.m
%
% Author: Michele Ducceschi, March 2026
%==========================================================================

clear; close all; clc;

%% ========================== USER PARAMETERS ==============================

% --- Paths ----------------------------------------------------------------
compMatPath = './polymaxRIR_results/comparison/comparison.mat';
outDir      = './polymaxRIR_results/comparison';

% --- Room and specimen geometry -------------------------------------------
%   V_room    : room volume [m³]
%   c_sound   : speed of sound [m/s]
%   S_spec    : one-sided area of the absorption specimen [m²]
%                (the area that faces the room — do NOT double-count)
V_room  = 277;    % [m³]
c_sound = 343;    % [m/s]
S_spec  = 10;    % [m²]

% --- Band resolution for output -------------------------------------------
%   'fullOctave' | 'thirdOctave'
%   Must match bandMode used in analyseResults.m so that dispEdges aligns.
bandMode = 'thirdOctave';

% --- Amplitude-weighting --------------------------------------------------
%   When true, per-mode α estimates are weighted by the mean residue
%   amplitude |r_p| = mean_IRs( √(a_p² + b_p²) ).  This down-weights
%   poorly-excited or weakly-resolved modes.
useAmplitudeWeighting = true;

% --- Output ---------------------------------------------------------------
verbose = true;

%% ========================== LOAD comparison.mat ==========================

fprintf('\n=== computeAbsorption ===\n');
if ~exist(compMatPath, 'file')
    error('computeAbsorption:noComp', ...
        'comparison.mat not found: %s\n  Run analyseResults.m first.', compMatPath);
end

comp = load(compMatPath);

emptyBasis    = comp.emptyBasis;
specimenBasis = comp.specimenBasis;   % same fBasis, only cBasis differs
emptyFits     = comp.emptyFits;
specimenFits  = comp.specimenFits;
dispEdges     = comp.dispEdges(:).';
dispCentres   = comp.dispCentres(:).';
fTransition   = comp.fTransition;
NDispBands    = numel(dispCentres);

emptyMeanT60      = comp.emptyMeanT60(:).';
emptyStdT60       = comp.emptyStdT60(:).';
specimenMeanT60   = comp.specimenMeanT60(:).';
specimenStdT60    = comp.specimenStdT60(:).';
emptyMeanAIT20    = comp.emptyMeanAIT20(:).';
specimenMeanAIT20 = comp.specimenMeanAIT20(:).';

% Load the original per-file data arrays (needed for per-IR WLS)
fprintf('Loading source MAT files for per-IR regression...\n');
emptyData    = loadConditionMats(comp.emptyDir,    verbose);
specimenData = loadConditionMats(comp.specimenDir, verbose);
NEmpty    = numel(emptyData);
NSpecimen = numel(specimenData);

fprintf('Loaded  : %s\n', compMatPath);
fprintf('Modes   : %d basis poles\n', emptyBasis.NBasis);
fprintf('Bands   : %d  (%s)  |  [%.1f – %.1f] Hz\n', ...
    NDispBands, bandMode, dispEdges(1), dispEdges(end));
fprintf('fTrans  : %.1f Hz\n', fTransition);
fprintf('V_room  : %.1f m³  |  S_spec : %.2f m²\n', V_room, S_spec);

%% =================== PER-IR WEIGHTED REGRESSION =========================
%
%   For each basis pole p we have:
%     - NEmpty  observations of c_E,i,p  with weight w_E,i,p = |r_E,i,p|
%     - NSpecimen observations of c_S,j,p with weight w_S,j,p = |r_S,j,p|
%
%   Model:  c_S,j,p − c_E,i,p  =  K_abs · α_p  +  ε_ij
%   where   K_abs = c_sound · S_spec / (8 V)
%
%   For each mode p the WLS solution reduces to:
%       Δc_p^WLS = c̄_S,p(weighted) − c̄_E,p(weighted)
%       α_p      = Δc_p^WLS / K_abs
%
%   Standard error of α_p comes from the propagated uncertainty in the
%   two weighted means:
%       Var(Δc) ≈ Var(c̄_S) + Var(c̄_E)
%       Var(c̄)  = Σ w_i²·(c_i − c̄)² / (Σ w_i)²
%
%   A spatial consistency index R_p measures how stable the amplitude
%   ratio |r_S,p|/|r_E,p| is across IRs.  R_p → 1 means the specimen
%   did not perturb mode shape; R_p ≪ 1 flags poor spatial coverage
%   (mode near a nodal line in one condition).  R_p is used as a
%   second-stage reliability weight on top of the amplitude weight.

K_abs = c_sound * S_spec / (8 * V_room);   % Δc = K_abs · α

% --- Collect per-IR damping and amplitudes for every basis pole -----------
%   cMat_E(p,i) = damping of pole p in empty IR i  (NaN if not found)
%   rMat_E(p,i) = |residue| of pole p in empty IR i
fBasis   = emptyBasis.fBasis(:);
cBasis_E = emptyBasis.cBasis(:);
cBasis_S = specimenBasis.cBasis(:);
NBasis   = emptyBasis.NBasis;
binWidth = emptyBasis.binWidth;

cMat_E = nan(NBasis, NEmpty);
rMat_E = nan(NBasis, NEmpty);
cMat_S = nan(NBasis, NSpecimen);
rMat_S = nan(NBasis, NSpecimen);

for iF = 1 : NEmpty
    fAll = emptyData(iF).fAll(:);
    cAll = emptyData(iF).cAll(:);
    aAll = emptyFits.aBasis{iF};
    bAll = emptyFits.bBasis{iF};
    for ip = 1 : NBasis
        [dMin, idx] = min(abs(fAll - fBasis(ip)));
        if dMin <= binWidth
            cMat_E(ip, iF) = cAll(idx);
        end
        rMat_E(ip, iF) = sqrt(aAll(ip)^2 + bAll(ip)^2);
    end
end

for iF = 1 : NSpecimen
    fAll = specimenData(iF).fAll(:);
    cAll = specimenData(iF).cAll(:);
    aAll = specimenFits.aBasis{iF};
    bAll = specimenFits.bBasis{iF};
    for ip = 1 : NBasis
        [dMin, idx] = min(abs(fAll - fBasis(ip)));
        if dMin <= binWidth
            cMat_S(ip, iF) = cAll(idx);
        end
        rMat_S(ip, iF) = sqrt(aAll(ip)^2 + bAll(ip)^2);
    end
end

% --- Per-mode ALL-PAIRS WLS regression ------------------------------------
%
%   For mode p we form all NEmpty × NSpecimen cross-pairs (i,j).
%   Each pair is one observation:
%       Δc_{ij} = c_{S,j,p} − c_{E,i,p}   with weight  w_{ij} = w_{E,i} · w_{S,j}
%
%   The single-parameter WLS solution for α_p:
%       α_p = (Σ_{ij} w_{ij} · Δc_{ij}) / (K_abs · Σ_{ij} w_{ij})
%
%   Residual-based standard error (honest: uses actual fit residuals):
%       SE(α_p)² = Σ_{ij} w_{ij}²·(Δc_{ij} − K_abs·α_p)²
%                  / (K_abs² · (Σ_{ij} w_{ij})²)
%
%   This uses O(N_E × N_S) equations per mode rather than O(N_E + N_S),
%   making full use of all source-receiver combinations.  The product weight
%   w_{ij} = w_{E,i} · w_{S,j} automatically down-weights pairs where either
%   IR fails to resolve the mode (near a nodal line or low SNR).
%
%   Outlier robustness: pairs whose |Δc_{ij} − K_abs·α_p| exceeds
%   OUTLIER_SIGMA · median_abs_deviation are flagged and excluded from the
%   SE estimate (they still enter the mean for stability).

OUTLIER_SIGMA = 3.0;

deltaC_wls   = nan(NBasis, 1);
deltaC_se    = nan(NBasis, 1);
alphaMod_wls = nan(NBasis, 1);
alphaMod_se  = nan(NBasis, 1);
consistencyR = nan(NBasis, 1);
rWeight      = nan(NBasis, 1);
nPairsUsed   = zeros(NBasis, 1);

for ip = 1 : NBasis
    wE = rMat_E(ip, :);  wE(isnan(wE)) = 0;
    wS = rMat_S(ip, :);  wS(isnan(wS)) = 0;
    cE = cMat_E(ip, :);
    cS = cMat_S(ip, :);

    okE = wE > 0 & ~isnan(cE);
    okS = wS > 0 & ~isnan(cS);
    if ~any(okE) || ~any(okS), continue; end

    wE_ = wE(okE)';   cE_ = cE(okE)';   % column vectors
    wS_ = wS(okS)';   cS_ = cS(okS)';

    nE = numel(cE_);  nS = numel(cS_);

    % All-pairs outer products: NE×NS matrices
    Wmat   = wE_ * wS_';          % pair weights
    DeltaC = cS_' - cE_;          % Δc_{ij} = c_S,j − c_E,i  (NE×NS)

    w_flat = Wmat(:);
    dc_flat = DeltaC(:);

    W_sum    = sum(w_flat);
    if W_sum < eps, continue; end

    % WLS estimate of α_p
    alphaP = sum(w_flat .* dc_flat) / (K_abs * W_sum);

    % Residuals and robust SE
    resid   = dc_flat - K_abs * alphaP;
    madR    = median(abs(resid - median(resid)));
    inliers = abs(resid) <= OUTLIER_SIGMA * (madR / 0.6745 + eps);
    w_in    = w_flat(inliers);
    r_in    = resid(inliers);
    W_in    = sum(w_in);
    if W_in > 0
        seP = sqrt(sum((w_in.^2) .* (r_in.^2))) / (K_abs * W_in);
    else
        seP = NaN;
    end

    deltaC_wls(ip)   = alphaP * K_abs;
    deltaC_se(ip)    = seP   * K_abs;
    alphaMod_wls(ip) = alphaP;
    alphaMod_se(ip)  = seP;
    nPairsUsed(ip)   = nE * nS;

    % Spatial consistency: coefficient of variation of per-IR |r|,
    % penalises modes that are inconsistently excited across positions.
    cvE = std(wE_) / (mean(wE_) + eps);
    cvS = std(wS_) / (mean(wS_) + eps);
    consistencyR(ip) = exp(-(cvE^2 + cvS^2));

    % Combined reliability: geometric mean amplitude × spatial consistency
    rWeight(ip) = sqrt(mean(wE_) * mean(wS_)) * consistencyR(ip);
end

rWeight = rWeight / (max(rWeight, [], 'omitnan') + eps);

% Hard-clamp α to [0,1]
alphaMod_raw_wls = alphaMod_wls;
alphaMod_wls     = max(0, alphaMod_wls);   % no upper clamp — values > 1 are physically informative

% Keep simple basis-median estimate for reference
rMag         = mean(rMat_E, 2, 'omitnan');
rMag         = rMag / (max(rMag) + eps);
deltaC       = cBasis_S - cBasis_E;
alphaMod     = max(0, deltaC / K_abs);
alphaMod_raw = deltaC / K_abs;

if verbose
    fprintf('\nAll-pairs regression: max pairs/mode = %d  (N_E=%d × N_S=%d)\n', ...
        NEmpty*NSpecimen, NEmpty, NSpecimen);
    nBelow = nnz(alphaMod_raw_wls < 0);
    nAbove = nnz(alphaMod_raw_wls > 1);
    if nBelow + nAbove > 0
        fprintf('Clamping: %d below 0 | %d above 1 (of %d modes)\n', ...
            nBelow, nAbove, NBasis);
    end
end

%% ============= LF: BAND-AGGREGATED α WITH UNCERTAINTY ===================
%
%   Aggregate per-mode WLS α estimates within each display band.
%   Weight by combined reliability (amplitude × spatial consistency).
%   std reflects mode-to-mode spread within the band.

alphaLF     = nan(1, NDispBands);
alphaLF_std = nan(1, NDispBands);
nModesLF    = zeros(1, NDispBands);

% Also compute the simpler median-based band average for comparison
alphaLF_med = nan(1, NDispBands);

nBelow = nnz(alphaMod_raw_wls < 0);
nAbove = nnz(alphaMod_raw_wls > 1);
if verbose && (nBelow + nAbove) > 0
    fprintf('\nPer-mode α clamping (WLS): %d below 0 | %d above 1 (out of %d)\n', ...
        nBelow, nAbove, NBasis);
end

for ib = 1 : NDispBands
    if dispCentres(ib) >= fTransition, continue; end
    inBand = (fBasis >= dispEdges(ib)) & (fBasis < dispEdges(ib + 1));
    if ~any(inBand), continue; end

    aVals = alphaMod_wls(inBand);
    wVals = rWeight(inBand);
    seVals = alphaMod_se(inBand) / K_abs;  % already in α units via /K_abs above

    ok = ~isnan(aVals) & ~isnan(wVals) & wVals > 0;
    if ~any(ok), continue; end

    aVals = aVals(ok);  wVals = wVals(ok);

    wNorm       = wVals / sum(wVals);
    alphaLF(ib) = sum(wNorm .* aVals);
    % Weighted std: mode-to-mode spread (dominant) + propagated SE
    alphaLF_std(ib) = sqrt(sum(wNorm .* (aVals - alphaLF(ib)).^2));
    nModesLF(ib) = nnz(ok);

    % Median-based for reference
    alphaLF_med(ib) = median(alphaMod(inBand), 'omitnan');
end

%% ============= HF: SABINE BAND α =========================================
%
%   For HF bands (f >= fTransition), use band-averaged T60 via Sabine:
%
%       ΔA    = 55.3 · V / c  · (1/T60_S − 1/T60_E)
%       α_HF  = ΔA / S_spec
%
%   Done twice: once with PolyMax T60 (emptyMeanT60/specimenMeanT60) and
%   once with AI T20 (emptyMeanAIT20/specimenMeanAIT20).  The PolyMax
%   estimate uses the modal damping directly; the AI T20 is the ISO-style
%   Schroeder estimate for cross-validation.

K = 55.3 * V_room / c_sound;   % Sabine constant [m³·s / m] = [m²·s]

% --- PolyMax-based ---
alphaHF_PM     = nan(1, NDispBands);
alphaHF_PM_std = nan(1, NDispBands);
for ib = 1 : NDispBands
    if dispCentres(ib) < fTransition, continue; end
    T60_E = emptyMeanT60(ib);
    T60_S = specimenMeanT60(ib);
    if isnan(T60_E) || isnan(T60_S) || T60_E <= 0 || T60_S <= 0, continue; end
    deltaA          = K * (1/T60_S - 1/T60_E);
    alphaHF_PM(ib)  = max(0, deltaA / S_spec);
    % Uncertainty from T60 std via Gaussian error propagation:
    % σ_α ≈ (K/S_spec) · √( (σ_S/T60_S²)² + (σ_E/T60_E²)² )
    sE = emptyStdT60(ib);   sS = specimenStdT60(ib);
    if ~isnan(sE) && ~isnan(sS)
        alphaHF_PM_std(ib) = (K / S_spec) * ...
            sqrt((sS / T60_S^2)^2 + (sE / T60_E^2)^2);
    end
end

% --- AI T20-based (all bands) ---
alphaAI     = nan(1, NDispBands);
for ib = 1 : NDispBands
    T60_E = emptyMeanAIT20(ib);
    T60_S = specimenMeanAIT20(ib);
    if isnan(T60_E) || isnan(T60_S) || T60_E <= 0 || T60_S <= 0, continue; end
    deltaA      = K * (1/T60_S - 1/T60_E);
    alphaAI(ib) = max(0, deltaA / S_spec);
end
alphaHF_AI = alphaAI;   % keep for backward-compat in saved struct

%% ============= COMBINE LF + HF INTO FINAL α CURVE =======================
%
%   Below fTransition : use per-mode weighted estimate (LF modal method)
%   Above fTransition : use Sabine PolyMax T60 (HF diffuse method)
%   Std is taken from whichever regime applies.

alphaFinal     = nan(1, NDispBands);
alphaFinal_std = nan(1, NDispBands);

for ib = 1 : NDispBands
    if dispCentres(ib) < fTransition
        alphaFinal(ib)     = alphaLF(ib);
        alphaFinal_std(ib) = alphaLF_std(ib);
    else
        alphaFinal(ib)     = alphaHF_PM(ib);
        alphaFinal_std(ib) = alphaHF_PM_std(ib);
    end
end

%% ========================= CONSOLE OUTPUT ================================

if verbose
    fprintf('\n%-8s  %8s  %8s  %8s  %8s  %7s  %s\n', ...
        'fc(Hz)', 'α_final', 'σ_α', 'α_HF_AI', 'α_HF_PM', 'N_modes', 'Regime');
    fprintf('%s\n', repmat('-', 1, 68));
    for ib = 1 : NDispBands
        if dispCentres(ib) < fTransition
            reg = 'LF';
        else
            reg = 'HF';
        end
        fprintf('%8.0f  %8.4f  %8.4f  %8.4f  %8.4f  %7d  %s\n', ...
            dispCentres(ib), ...
            alphaFinal(ib),     alphaFinal_std(ib), ...
            alphaAI(ib),        alphaHF_PM(ib), ...
            nModesLF(ib),       reg);
    end
end

%% ============================= FIGURES ===================================

colLF   = [0.18 0.40 0.75];   % blue  — LF modal WLS
colLFm  = [0.55 0.70 0.95];   % light blue — LF modal median (reference)
colHF   = [0.82 0.20 0.15];   % red   — HF Sabine
colAI   = [0.50 0.50 0.50];   % grey  — AI T20 Sabine (ISO)
colSch  = [0.10 0.55 0.10];   % green — Schroeder transition

%--------------------------------------------------------------------------
%  Figure 1 — Final α curve: all methods overlaid
%--------------------------------------------------------------------------
hFig = figure('Name', 'Absorption coefficient', ...
    'Color', 'w', 'Position', [100 100 1050 560]);
hold on;

yLim_guess = max(1.15, max([alphaLF, alphaHF_PM, alphaAI], [], 'omitnan') * 1.1 + 0.05);
if isnan(yLim_guess) || yLim_guess <= 0, yLim_guess = 1.15; end
patch([dispEdges(1), fTransition, fTransition, dispEdges(1)], ...
    [0, 0, yLim_guess, yLim_guess], ...
    [0.92 0.96 0.92], 'EdgeColor', 'none', 'FaceAlpha', 0.5, ...
    'HandleVisibility', 'off');

% AI T20 Sabine — full range (ISO reference)
hasAI = ~isnan(alphaAI);
if any(hasAI)
    plot(dispCentres(hasAI), alphaAI(hasAI), 's--', ...
        'Color', colAI, 'MarkerFaceColor', colAI, ...
        'MarkerSize', 7, 'LineWidth', 1.2, ...
        'DisplayName', 'Sabine / AI T_{20}  (ISO reference)');
end

% LF: median-based (light blue, thin) — simple reference
hasLFm = ~isnan(alphaLF_med) & (dispCentres < fTransition);
if any(hasLFm)
    plot(dispCentres(hasLFm), alphaLF_med(hasLFm), 'o:', ...
        'Color', colLFm, 'MarkerFaceColor', colLFm, ...
        'MarkerSize', 6, 'LineWidth', 1.0, ...
        'DisplayName', 'LF — median \Deltac  (basis only)');
end

% HF: PolyMax Sabine
hasHF = ~isnan(alphaHF_PM) & (dispCentres >= fTransition);
if any(hasHF)
    errorbar(dispCentres(hasHF), alphaHF_PM(hasHF), alphaHF_PM_std(hasHF), ...
        's-', 'Color', colHF, 'MarkerFaceColor', colHF, ...
        'MarkerSize', 7, 'LineWidth', 1.8, 'CapSize', 6, ...
        'DisplayName', 'HF — Sabine / PolyMax T_{60}  (mean \pm 1\sigma)');
end

% LF: WLS (blue, main result)
hasLF = ~isnan(alphaLF) & (dispCentres < fTransition);
if any(hasLF)
    errorbar(dispCentres(hasLF), alphaLF(hasLF), alphaLF_std(hasLF), ...
        'o-', 'Color', colLF, 'MarkerFaceColor', colLF, ...
        'MarkerSize', 7, 'LineWidth', 1.8, 'CapSize', 6, ...
        'DisplayName', 'LF — per-IR WLS \Deltac  (mean \pm 1\sigma)');
    for ib = find(hasLF)
        if nModesLF(ib) > 0
            text(dispCentres(ib), alphaLF(ib) + alphaLF_std(ib) + 0.025, ...
                sprintf('%d', nModesLF(ib)), ...
                'HorizontalAlignment', 'center', 'FontSize', 7, 'Color', colLF);
        end
    end
end

yline(1, 'k--', 'LineWidth', 1.2, 'Label', '\alpha = 1', ...
    'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
xline(fTransition, '-', 'Color', colSch, 'LineWidth', 1.5, ...
    'Label', 'f_{Sch}', 'LabelVerticalAlignment', 'top', 'HandleVisibility', 'off');
text(fTransition * 0.55, min(yLim_guess * 0.96, 1.08), 'LF  (modal \Deltac)', ...
    'HorizontalAlignment', 'right', 'Color', [0.3 0.55 0.3], 'FontSize', 9);
text(fTransition * 1.45, min(yLim_guess * 0.96, 1.08), 'HF  (Sabine)', ...
    'HorizontalAlignment', 'left',  'Color', [0.3 0.55 0.3], 'FontSize', 9);

set(gca, 'XScale', 'log');
xlabel('Band centre frequency  f_c  (Hz)', 'FontSize', 12);
ylabel('Sound absorption coefficient  \alpha', 'FontSize', 12);
xlim([dispEdges(1), dispEdges(end)]);  ylim([0, yLim_guess]);
legend('show', 'Location', 'best', 'FontSize', 10);
grid on;  box on;
title({'Sound absorption coefficient  \alpha_{spec}', ...
    sprintf('V = %.0f m³  |  S_{spec} = %.2f m²  |  %s bands', ...
    V_room, S_spec, bandMode)}, 'FontSize', 12);

%--------------------------------------------------------------------------
%  Figure 2 — Per-mode LF scatter: α vs f, coloured by reliability
%--------------------------------------------------------------------------
lfMask = fBasis < fTransition & ~isnan(alphaMod_wls) & ~isnan(rWeight);

hFig2 = figure('Name', 'Per-mode absorption', ...
    'Color', 'w', 'Position', [100 720 1100 480]);

subplot(1, 2, 1);
hold on;
sc = scatter(fBasis(lfMask), alphaMod_wls(lfMask), 30, ...
    rWeight(lfMask), 'filled', 'MarkerFaceAlpha', 0.75);
colormap(gca, parula);  cb = colorbar;
cb.Label.String = 'Reliability  (amplitude \times consistency)';
clim([0, 1]);
yline(1, 'k--', 'LineWidth', 1.2, 'Label', '\alpha = 1', ...
    'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
yline(0, 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');
xline(fTransition, '--', 'Color', colSch, 'LineWidth', 1.2, ...
    'Label', 'f_{Sch}', 'HandleVisibility', 'off');
xlabel('Modal frequency  f_p  (Hz)');
ylabel('\alpha_{spec,p}  (per mode)');
title({'Per-mode \alpha_{spec}  (all-pairs WLS)', ...
    sprintf('up to %d pairs/mode  (N_E=%d × N_S=%d)', ...
    NEmpty*NSpecimen, NEmpty, NSpecimen)});
grid on;  box on;

subplot(1, 2, 2);
hold on;
% Reliability histogram
histogram(rWeight(lfMask), 25, 'FaceColor', colLF, 'EdgeColor', 'none', ...
    'Normalization', 'probability');
xline(0.5, 'r--', 'LineWidth', 1.2, ...
    'Label', 'reliability = 0.5', 'HandleVisibility', 'off');
xlabel('Reliability index');
ylabel('Probability');
title({'Reliability distribution', sprintf('%d LF modes below f_{Sch}', nnz(lfMask))});
grid on;  box on;

%% ============================= SAVE ======================================

if ~exist(outDir, 'dir'), mkdir(outDir); end

savefig(hFig, fullfile(outDir, 'absorption_alpha.fig'));
fprintf('\nFigure saved : %s\n', fullfile(outDir, 'absorption_alpha.fig'));

absS.alphaFinal       = alphaFinal;
absS.alphaFinal_std   = alphaFinal_std;
absS.alphaLF          = alphaLF;
absS.alphaLF_std      = alphaLF_std;
absS.alphaHF_PM       = alphaHF_PM;
absS.alphaHF_PM_std   = alphaHF_PM_std;
absS.alphaAI          = alphaAI;          % Sabine/AI T20, all bands
absS.alphaMod         = alphaMod;         % per-mode α (LF only)
absS.alphaMod_raw     = alphaMod_raw;     % unclamped per-mode α
absS.fBasis           = fBasis;
absS.deltaC           = deltaC;
absS.rMag             = rMag;
absS.nPairsUsed       = nPairsUsed;       % NE×NS pairs used per mode
absS.dispEdges        = dispEdges;
absS.dispCentres      = dispCentres;
absS.fTransition      = fTransition;
absS.V_room           = V_room;
absS.c_sound          = c_sound;
absS.S_spec           = S_spec;
absS.bandMode         = bandMode;
absS.useAmplitudeWeighting = useAmplitudeWeighting;
absS.processedDate    = datestr(now);  %#ok<TNOW1,DATST>

savePath = fullfile(outDir, 'absorption.mat');
save(savePath, '-struct', 'absS');
fprintf('MAT  saved  : %s\n', savePath);

fprintf('\n=== DONE ===\n');

%% ========================= LOCAL FUNCTION ================================

function data = loadConditionMats(folder, verbose)
%LOADCONDITIONMATS  Load all *_polymaxRIR.mat files from FOLDER.
matFiles = dir(fullfile(folder, '*_polymaxRIR.mat'));
if isempty(matFiles)
    error('computeAbsorption:noMats', ...
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

% -------------------------------------------------------------------------

function s = conditional(cond, a, b)
if cond, s = a; else, s = b; end
end