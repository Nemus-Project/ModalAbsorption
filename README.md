# Modal Absorption Coefficients

Modal analysis, synthesis, and absorption measurement of room impulse responses using spectral peak-picking with independent half-power bandwidth damping.

**Author:** Michele Ducceschi, March 2026

---

## Overview

`modalRIR` is a MATLAB pipeline for physically-grounded analysis of room impulse responses (RIRs). It estimates modal parameters — frequencies **f_p**, damping coefficients **c_p**, and complex residues **r_p = a_p + j b_p** — and reconstructs a causal synthetic RIR as a direct modal sum:

```
h_syn(n) = 2 Re { Σ_p  r_p  exp(λ_p · n/fs) },    λ_p = −c_p + j√(ω_p² − c_p²)
```

Synthesis is done directly in the time domain — not via IFFT — so decay tails are physically correct regardless of RIR length.

The analysis range is split at the **Schroeder frequency** f_Sch into two physically motivated regimes:

| Regime | Range | Method | Damping estimate |
|--------|-------|--------|-----------------|
| **LF** | f < f_Sch | Spectral peak-picking on \|H(f)\| | Half-power (−3 dB) bandwidth per mode |
| **HF** | f ≥ f_Sch | Spectral peak-picking on \|H(f)\| | Uniform Schroeder T30 per band |

The LF damping estimate is completely independent of the Schroeder integral, T20, or T30. When the −3 dB crossing is blocked by a neighbouring peak or band edge, a local Lorentzian least-squares fit is used instead.

f_Sch is computed from the mid-frequency T20 (500–1000 Hz):

```
f_Sch = 2000 √(T20_mid / V_room)
```

---

## Files

| File | Role |
|------|------|
| `modalRIR.m` | Main batch script — modal analysis, synthesis, per-file `.mat` and `.wav` output |
| `acoustic_indices.m` | Room-acoustic indices (T20, C80, Ts, BR, TR) — standalone reusable function |
| `analyseResults.m` | Post-processing comparison of two conditions (e.g. empty vs specimen) |
| `computeAbsorption.m` | Sound absorption coefficient estimation from the comparison results |

All pipeline helper functions are local functions inside their respective scripts.

---

## Pipeline

```
measured_IRs/
  empty/*.wav          ──► modalRIR.m ──► modalRIR_results/empty/*_modalRIR.mat
  specimen/*.wav       ──► modalRIR.m ──► modalRIR_results/specimen/*_modalRIR.mat
                                                        │
                                                 analyseResults.m
                                                        │
                                           comparison/comparison.mat
                                           comparison/fig1–fig8.fig
                                           comparison/empty/*_commonBasis.mat
                                           comparison/empty/*_commonBasis_SYNTH.wav
                                           comparison/specimen/*_commonBasis.mat
                                           comparison/specimen/*_commonBasis_SYNTH.wav
                                                        │
                                            computeAbsorption.m
                                                        │
                                           comparison/absorption.mat
                                           comparison/absorption_alpha.fig
```

---

## Dependencies

| MATLAB Toolbox | Functions used |
|----------------|---------------|
| Signal Processing Toolbox | `filtfilt`, `designfilt`, `findpeaks` |
| Audio Toolbox | `audioread`, `audiowrite` |

---

## Step 1 — modalRIR.m

### Quick start

1. Place `.wav` files in a folder (one condition per folder).
2. Edit the **USER PARAMETERS** block at the top of `modalRIR.m`.
3. Run the script. All results are saved to `outRootDir/<folderName>/`.

### Batch mode

`modalRIR` processes all `.wav` files found in `inputFolder` in a single run. Each file is processed independently inside a `try/catch` block. A batch summary table is printed and saved as `batchSummary.mat`.

```matlab
inputFolder = './measured_IRs/empty';
outRootDir  = './modalRIR_results/';
```

### Two-band-grid architecture

**Processing bands** are narrow linear bands of fixed width `procBandWidth_Hz` [Hz]. Peak-picking runs on these. A constant bandwidth keeps the number of peaks per band manageable across the whole frequency range.

**Display bands** are ISO 266 / IEC 61260 octave or third-octave bands. All statistics, figures, and the saved `poleStats` struct use these. Poles are assigned to display bands by frequency after estimation.

Both grids are independent: `bandMode` and `procBandWidth_Hz` can be tuned separately.

### User parameters

#### File and output

```matlab
inputFolder = './measured_IRs/St Gobain';   % folder containing .wav files
outRootDir  = './modalRIR_results/';       % root for output subfolders
```

#### Analysis frequency range

```matlab
fLow  = 20;      % [Hz]  lower analysis limit
fHigh = 12000;   % [Hz]  upper analysis limit
```

#### Display bands

```matlab
bandMode        = 'thirdOctave';    % 'thirdOctave' | 'fullOctave' | 'custom'
customBandEdges = [50 100 200 ...]; % only used when bandMode = 'custom'
```

#### Processing bands

```matlab
procBandWidth_Hz = 35;   % [Hz]  linear processing band width
```

Typical range: 20–100 Hz.

#### Room geometry

```matlab
V_room  = 277;   % [m³]  room volume
c_sound = 343;   % [m/s] speed of sound
```

Used for Weyl mode density and Schroeder frequency estimation.

#### Transition frequency

```matlab
fTransition = [];   % [Hz]  set to [] for automatic
```

When empty, computed from the mid-frequency T20 from `acoustic_indices`. Set a numeric value to override.

#### LF peak-picking parameters

```matlab
LFpeakProm_dB  = 4.0;    % [dB]  minimum peak prominence
LFpeakDist_Hz  = 0.25;   % [Hz]  minimum inter-peak spacing
```

Controls which spectral peaks are retained as modes below f_Sch. Prominence is evaluated relative to the local spectral floor.

#### HF peak-picking parameters

```matlab
minPeakProm_dB = 0.5;   % [dB]  minimum peak prominence
minPeakDist_Hz = 0.5;   % [Hz]  minimum inter-peak spacing
```

#### Gain outlier rejection

```matlab
enableGainOutlierRejection = true;
gainOutlierMethod          = 'mad';   % 'mad' | 'std'
gainOutlierZthresh         = 4.0;
maxGainOutlierLSIters      = 5;
minModesToKeepAfterReject  = 2;
```

After the LS amplitude solve, poles with residue amplitudes more than `gainOutlierZthresh` robust Z-scores from the band median are discarded and the solve is repeated.

#### Iterative refinement

```matlab
useIterativeRefinement = true;
nRefinementIters       = 8;
```

When enabled, re-estimates poles from the spectral residual of the previous iteration, adding modes that the initial pass missed.

#### FIR correction

```matlab
useFIRCorrection = false;
```

When enabled, a linear-phase FIR filter is fitted to the residual between reference and synthetic signal to correct for broadband spectral mismatch.

#### Synthesis quality threshold

```matlab
nRMSE_time_keep_threshold = 5e-3;
```

Files with time-domain nRMSE above this value are flagged in the batch summary (but still saved).

### Output files (per IR)

| File | Description |
|------|-------------|
| `<n>_modalRIR.mat` | Full results struct (see fields below) |
| `<n>_modalRIR_SYNTH.wav` | Direct modal synthesis (peak-normalised) |
| `<n>_spectrum.fig` | Spectrum + residual + time-domain comparison |
| `<n>_bandStats.fig` | Per-band mode count, T60, and spacing |
| `<n>_LFcomparison.fig` | Modal T60 vs acoustic-indices T20 for LF bands |

### MAT file struct fields

| Field | Description |
|-------|-------------|
| `fAll`, `cAll`, `aAll`, `bAll` | All poles: frequency [Hz], damping [rad/s], residue components |
| `bandIdxAll` | Display-band index for each pole |
| `poleStats` | Per-display-band statistics struct |
| `globalStats` | Wideband summary struct |
| `refAI` | Full `acoustic_indices` output struct |
| `fTransition` | Schroeder frequency used [Hz] |
| `dispEdgesVec` | Display band edges [Hz] |
| `nRMSE_IIR`, `nRMSE_final` | Time-domain normalised RMSE |
| `fs` | Sample rate [Hz] |

### poleStats struct fields (per display band)

| Field | Description |
|-------|-------------|
| `flow`, `fhigh`, `fcBand` | Band edges and ISO centre [Hz] |
| `N` | Number of poles in band |
| `regime` | `'LF'` or `'HF'` |
| `modalDensity_Hz` | Measured modal density [modes/Hz] |
| `densityRatio` | Measured / Weyl theoretical density |
| `meanT60`, `stdT60`, `cvT60` | T60 statistics [s] |
| `meanC`, `stdC` | Damping coefficient statistics [rad/s] |
| `meanF`, `stdF` | Modal frequency statistics [Hz] |
| `meanSpacing_Hz`, `cvSpacing` | Inter-mode spacing statistics [Hz] |

---

## Step 2 — analyseResults.m

Compares modal parameters across two conditions (empty room vs room with absorption specimen). Loads all `*_modalRIR.mat` files from two result folders.

### Cache-skip behaviour

On each run, `analyseResults` checks whether `comparison.mat` and all per-IR `_commonBasis.mat` files already exist in `outDir`, and whether the cached MAT contains all required fields (including newer ones like C80 and mode-count matrices). If anything is missing or outdated, the full pipeline re-runs automatically. Delete `comparison.mat` to force a full re-run.

### Algorithm

1. **Common pole basis**: vote on shared modal frequencies from the empty-room IRs. A frequency bin enters the basis when at least `minFileFraction` of IRs contribute a pole. Representative `(f, c)` per basis pole = median across contributing IRs.
2. **Specimen damping**: for each basis frequency, find the nearest pole in each specimen IR within `binWidth_Hz` and take the median `c` across IRs. Frequencies are NOT re-voted — they are geometry-locked.
3. **Amplitude re-fit**: banded linear LS solve for residues `(a_p, b_p)` on the shared `(f, c)` basis, per IR.
4. **Statistics**: per-display-band T60, C80, mode count, modal density, and scalar metrics (BR, TR, Ts) from both the original per-IR pole sets and the common basis.

### User parameters

```matlab
emptyDir    = './modalRIR_results/empty';
specimenDir = './modalRIR_results/specimen';
outDir      = './modalRIR_results/comparison';

binWidth_Hz     = 0.5;    % [Hz]  frequency bin width for pole voting
minFileFraction = 0.4;    % majority vote threshold [0, 1]
bandMode        = 'thirdOctave';
lsBandWidth_Hz  = 50;     % [Hz]  linear LS band width for amplitude re-fit
overlapFrac     = 0.10;   % LS window extension beyond band edges
MinSecondsForFFT = 10;    % [s]   FFT zero-padding floor
```

### Figures produced

| Figure | Contents |
|--------|----------|
| `fig1_T60distributions.fig` | T60 histograms (LF and HF), empty vs specimen |
| `fig2_T60spaghetti.fig` | Per-IR T60 curves (faint) + mean (solid), both conditions in one panel |
| `fig3_T60comparison.fig` | Mean ± 1σ T60 per band: modal (solid) and AI T20 (dashed), both conditions |
| `fig4_spaghettiOverlay.fig` | Per-IR spaghetti with final mean ± 1σ overlay (modal + AI T20), one panel per condition |
| `fig5_heatmapGrid.fig` | Multi-metric heatmap: IRs × bands × {T60, C80, mode count} |
| `fig6_radarConsistency.fig` | Spatial consistency radar: per-IR T60, C80, BR, TR, Ts |
| `fig7_correlationMatrix.fig` | Pairwise correlation scatter matrix: T60, C80, modal density |
| `fig8_stripPlot.fig` | Per-IR strip plot: jittered dots for T60, C80, mode count |

### Output files

| File | Description |
|------|-------------|
| `comparison.mat` | Full comparison struct (bases, fits, T60/C80/mode-count matrices, scalar metrics, parameters) |
| `empty/<n>_commonBasis.mat` | Per-IR LS amplitudes on shared basis |
| `empty/<n>_commonBasis_SYNTH.wav` | Resynthesised IR on common basis (peak-normalised) |
| `specimen/<n>_commonBasis.mat` | As above for specimen condition |
| `specimen/<n>_commonBasis_SYNTH.wav` | As above for specimen condition |

---

## Step 3 — computeAbsorption.m

Estimates the sound absorption coefficient α of the specimen per octave or third-octave band from the `comparison.mat` output of `analyseResults`.

### Method

**LF regime (f < f_Sch) — per-mode all-pairs WLS regression**

Every (empty IR i, specimen IR j) pair is one observation of the modal damping shift:

```
Δc_{ij,p} = c_{S,j,p} − c_{E,i,p}   →   K_abs · α_p
```

with pair weight `w_{ij} = w_{E,i} · w_{S,j}` where `w = |r| = √(a² + b²)`. This all-pairs formulation uses N_E × N_S equations per mode, maximising statistical efficiency.

```
α_p = Σ_{ij} w_{ij} · Δc_{ij} / (K_abs · Σ_{ij} w_{ij})

K_abs = c_sound · S_spec / (8 V)
```

Per-band α is the reliability-weighted mean of α_p over all modes in the band, where the reliability weight combines amplitude magnitude and spatial consistency.

**HF regime (f ≥ f_Sch) — Sabine formula on modal T60**

```
ΔA    = 55.3 V / c · (1/T60_S − 1/T60_E)
α_HF  = ΔA / S_spec
```

Uncertainty is propagated analytically from the inter-IR T60 standard deviations.

**Advantages over ISO Schroeder / Sabine**

| Issue | ISO method | This method |
|-------|-----------|-------------|
| Nonlinearity bias | Averages T60 then inverts | Works directly in c (linear) |
| Noise floor | Schroeder tail corrupts late decay | Peak-picking robust to late-decay noise |
| Room absorption | Must measure or subtract | Cancels exactly via Δc |
| Uncertainty | Single T60 value | Per-mode spread → natural σ |
| Below f_Sch | T60 unreliable (few modes) | Mode-by-mode resolution |
| Spatial sampling | All IRs averaged equally | Amplitude-weighted by \|r_{ij}\| |

All α values are hard-clamped to [0, 1].

### User parameters

```matlab
compMatPath = './modalRIR_results/comparison/comparison.mat';
outDir      = './modalRIR_results/comparison';

V_room   = 277;    % [m³]  room volume
c_sound  = 343;    % [m/s] speed of sound
S_spec   = 10;     % [m²]  specimen area (one-sided, facing the room)
bandMode = 'thirdOctave';

useAmplitudeWeighting = true;   % weight per-mode α by |r_p|
```

### Figures produced

| Figure | Contents |
|--------|----------|
| `absorption_alpha.fig` | α per band: LF WLS (blue, errorbars), LF median (light blue), HF Sabine (red, errorbars), ISO reference Sabine/AI T20 (grey, full range) |
| *Per-mode absorption* | Scatter of α_p vs f_p coloured by reliability; reliability distribution histogram |

### Output files

| File | Description |
|------|-------------|
| `absorption.mat` | α curve, uncertainties, per-mode estimates, reliability weights, all parameters |

---

## Standalone usage

### acoustic_indices.m

```matlab
[h, fs] = audioread('myRIR.wav');
res = acoustic_indices(h, fs, bandsPerOctave, freqLims, tEarly, alignThreshold_dB);
```

| Input | Description |
|-------|-------------|
| `h` | Mono RIR vector |
| `fs` | Sample rate [Hz] |
| `bandsPerOctave` | 1 for octave, 3 for third-octave |
| `freqLims` | `[fMin fMax]` [Hz] |
| `tEarly` | Early/late boundary for C80 [s] (ISO 3382: 0.08 s) |
| `alignThreshold_dB` | Onset detection threshold below peak [dB] (e.g. −5) |

| Output field | Description |
|-------------|-------------|
| `bandCentersHz` | ISO band centre frequencies [Hz] |
| `T20` | Reverberation time per band [s] — ISO 3382 T20 estimator (fit −5 to −25 dB, extrapolated to 60 dB) |
| `C80` | Clarity per band [dB] |
| `Ts` | Wideband centre time [s] |
| `BR` | Bass ratio: (T20₁₂₅ + T20₂₅₀) / (T20₅₀₀ + T20₁₀₀₀) |
| `TR` | Treble ratio: (T20₂₀₀₀ + T20₄₀₀₀) / (T20₅₀₀ + T20₁₀₀₀) |
| `hAligned` | Onset-aligned, mid-T60-cropped RIR (column vector) |

---

## Algorithm notes

### Direct time-domain synthesis

Synthesis uses the modal sum directly rather than IFFT:

```
h_syn(n) = 2 Re { Σ_p  (a_p + j b_p)  exp(λ_p · n/fs) }
```

An IFFT-based approach treats the spectrum as periodic with period T = N/fs. Any mode whose T60 exceeds T wraps its tail back onto the start of the block, producing a spurious burst at the end. Direct synthesis has no such window constraint — each mode decays correctly to zero. Poles are processed in memory-safe chunks (~50 MB per chunk).

### LF damping: half-power bandwidth

For each detected spectral peak, the damping coefficient σ is estimated from the −3 dB bandwidth Δf₃dB:

```
σ = π · Δf₃dB     [rad/s]
T60 = 3·ln(10) / σ
```

Sub-bin interpolation is used for accuracy. When a neighbouring peak blocks the −3 dB crossing, a local Lorentzian LS fit is used instead:

```
|H(f)|² ≈ A / [ (f − f₀)² + γ² ]     →     σ = 2π · γ
```

This is an independent per-mode measurement that does not depend on Schroeder integration, T20, or T30.

### T20 as the T60 estimator

`acoustic_indices` uses exclusively the **T20** estimator throughout: the Schroeder decay curve is fitted over −5 to −25 dB and the slope is extrapolated to a 60 dB drop. Per ISO 3382, T20 is an unbiased estimate of T60 — no factor of 2 is applied. T20 is used instead of T30 because it requires only 25 dB of usable dynamic range, making it more robust for low-SNR or short IRs.

### Amplitude solve

Residues are found by banded LS from the partial-fraction model:

```
H(jω) = Σ_p  r_p / (jω − λ_p)  +  r_p* / (jω − λ_p*)
```

Solved per processing band over a window extended by `overlapFrac` to reduce boundary artefacts. Tikhonov regularisation with `ε = 10⁻⁸ × mean(diag(R'R))`.

### ISO band edge construction

```
fc(n) = 1000 · 2^(n/B)   Hz,    B = bandsPerOctave
lo(n) = fc(n) / 2^(1/2B),    hi(n) = fc(n) · 2^(1/2B)
```

Edges are always exactly N+1 values; `fLow` and `fHigh` select which bands to include without modifying their widths.

---

## References

Schroeder M. R. (1965).
*New method of measuring reverberation time.*
Journal of the Acoustical Society of America, **37**(3), 409–412.

ISO 3382-1:2009. *Acoustics — Measurement of room acoustic parameters — Part 1: Performance spaces.*

ISO 266:1997. *Acoustics — Preferred frequencies.*

IEC 61260-1:2014. *Electroacoustics — Octave-band and fractional-octave-band filters.*
