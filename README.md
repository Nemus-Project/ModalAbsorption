# polymaxRIR

Modal analysis and synthesis of room impulse responses using the PolyMax frequency-domain estimator.

**Author:** Michele Ducceschi, March 2026

---

## What it does

Given a measured room impulse response (RIR) as a `.wav` file, `polymaxRIR` estimates the room's modal parameters — frequencies, damping coefficients, and complex amplitudes — and reconstructs a synthetic RIR as a sum of exponentially decaying sinusoids.

The analysis range is split at the **Schroeder frequency** into two physically motivated regimes:

| Regime | Method | Rationale |
|--------|--------|-----------|
| **LF** `f < f_Sch` | PolyMax stability-diagram pole fitting | Modes are individually resolvable; precise estimation is possible |
| **HF** `f ≥ f_Sch` | Spectral peak-picking on \|H(f)\| | Dense modal overlap makes pole-by-pole resolution unreliable; peaks with Schroeder-T30 damping are the best available estimate |

---

## Files

| File | Role |
|------|------|
| `polymaxRIR.m` | Main script — full analysis and synthesis pipeline |
| `polymax.m` | PolyMax estimator — standalone reusable function |
| `acoustic_indices.m` | Room-acoustic indices — standalone reusable function |

All pipeline helper functions (band construction, parameter selection, LS amplitude solve, outlier rejection, statistics, plotting) live as local functions at the bottom of `polymaxRIR.m`.

---

## Dependencies

| MATLAB Toolbox | Used for |
|----------------|----------|
| Signal Processing Toolbox | `filtfilt`, `designfilt`, `findpeaks` |
| Statistics and Machine Learning Toolbox | `linkage`, `cluster` |
| Audio Toolbox | `audioread`, `audiowrite` |

---

## Quick start

1. Place your RIR as a `.wav` file somewhere accessible.
2. Open `polymaxRIR.m` and fill in the **USER PARAMETERS** block at the top.
3. Run the script. Results are saved to `outRootDir` (`./polymaxRIR_results/` by default).

---

## User parameters

### File

```matlab
filename  = 'Empty_S1R01';   % .wav filename without extension
inputPath = './';             % folder containing the file
```

### Analysis range and band mode

```matlab
fLow  = 20;      % [Hz]  lower analysis limit
fHigh = 4000;    % [Hz]  upper analysis limit

bandMode = 'thirdOctave';   % 'thirdOctave' | 'fullOctave' | 'custom'
```

Bands follow the exact ISO 266 / IEC 61260 centre frequencies. `fLow` and `fHigh` select which standard bands to include; the band edges themselves are never adjusted to match these values, so every band retains its standard ISO width.

For `'custom'` mode, provide explicit edges:

```matlab
customBandEdges = [50 100 200 400 800 1600 3200 4000];
```

### Room geometry

```matlab
V_room  = 200;   % [m³]  room volume
c_sound = 343;   % [m/s] speed of sound
```

Used to compute the Weyl modal density estimate (modes/Hz) that sizes the PolyMax polynomial order for each LF band.

### Transition (Schroeder) frequency

```matlab
fTransition = [];   % [Hz]  set to [] for automatic computation
```

When empty, `fTransition` is computed as:

```
f_Sch = 2000 * sqrt(T60_mid / V_room)
```

where `T60_mid` is the mean T60 in the 500–1000 Hz range as measured by `acoustic_indices`. Supply a numeric value to fix the transition manually.

### PolyMax parameters

```matlab
maxPolyOrderLF = 2000;   % hard ceiling on polynomial order (LF bands)
maxPolyIters   = 40;     % target number of order steps per band
```

The per-band order is set to `min(6 × N_expected, Nbins/4, maxPolyOrderLF, 900)`, where `N_expected` is the Weyl-theory mode count. The 900-sample absolute cap prevents intractably large companion matrices for wide octave bands at mid-frequencies.

```matlab
tolf_init = 0.10;   % relative frequency stability tolerance
told_init = 0.10;   % relative damping  stability tolerance
```

A pole is accepted as stable when both its relative frequency shift and relative damping shift between consecutive polynomial orders are below these thresholds.

### HF peak-picking

```matlab
minPeakProm_dB = 1.0;   % [dB]  minimum peak prominence above local floor
minPeakDist_Hz = 0.5;   % [Hz]  minimum distance between peaks
```

Above `fTransition`, peaks are detected in `|H(f)|` using `findpeaks`. Damping is set uniformly to `3·ln(10) / T60est` from the per-band Schroeder T30.

### T60 rejection ceiling

```matlab
T60cutMult  = 3.0;   % T60cut = T60cutMult × T60est  (per band)
T60cutFloor = 5.0;   % [s]  minimum T60cut regardless of T60est
```

Poles with `T60 = 3·ln(10)/c > T60cut` are discarded. This removes spurious low-damping poles that arise from noise or numerical artefacts.

### Gain outlier rejection

```matlab
enableGainOutlierRejection = true;
gainOutlierMethod          = 'mad';   % 'mad' (robust) | 'std' (classical)
gainOutlierZthresh         = 6.0;
```

After the LS amplitude solve, poles whose residue amplitudes deviate by more than `gainOutlierZthresh` robust Z-scores from the band median are discarded and amplitudes are re-solved. Iterates up to `maxGainOutlierLSIters` times.

### Iterative residual refinement *(optional)*

```matlab
useIterativeRefinement = false;
nRefinementIters       = 8;
```

When enabled, PolyMax is re-run on the spectral residual `H(f) - H_syn(f)` to find modes missed in the first pass. Useful for dense LF fields but substantially increases runtime. Stops early if nRMSE improvement per iteration falls below `minNRMSEImprovement`.

### FIR correction filter *(optional)*

```matlab
useFIRCorrection = false;
FIR_Ntaps        = 256;
```

A length-`FIR_Ntaps` least-squares FIR filter is fitted to compensate for systematic spectral residuals that the modal IIR model cannot represent (e.g. diffuse-field colouration). Applied as a post-processing step; the IIR output is unchanged.

### Acoustic indices

```matlab
bandsPerOctave    = 3;        % used for acoustic_indices analysis
freqLimsAI        = [20, fHigh];
tEarly            = 0.080;    % [s]  C80 early/late boundary
alignThreshold_dB = -5;       % [dB]  onset detection level below peak
```

`acoustic_indices` runs first on the raw RIR. Its outputs drive three things in the main pipeline: onset alignment, per-band T60 estimates for the T60cut and HF damping, and the automatic Schroeder frequency.

### Output

```matlab
outRootDir                = './polymaxRIR_results/';
nRMSE_time_keep_threshold = 5e-3;
verbose                   = true;
```

Files where `nRMSE_final > nRMSE_time_keep_threshold` are flagged `DISCARDED` in the console. No data is deleted; the flag is informational.

---

## Outputs

### Files written to `outRootDir`

| File | Contents |
|------|----------|
| `<filename>_polymaxRIR.mat` | All results (see table below) |
| `<filename>_polymaxRIR_SYNTH.wav` | Peak-normalised synthetic RIR |

### MAT struct fields

| Field | Description |
|-------|-------------|
| `fAll` | Modal frequencies [Hz], sorted ascending |
| `cAll` | Modal damping coefficients [rad/s];  T60 = 3·ln(10)/c |
| `aAll`, `bAll` | Real and imaginary residue components |
| `bandIdxAll` | Analysis band index each pole belongs to |
| `hFIR` | FIR correction filter taps (scalar 1 if unused) |
| `poleStats` | Per-band statistics struct array |
| `globalStats` | Global summary (total poles, mean T60, Schroeder freq) |
| `refAI` | Full `acoustic_indices` output struct |
| `nRMSE_IIR` | Normalised time-domain RMSE of the IIR-only synthesis |
| `nRMSE_final` | Normalised time-domain RMSE after FIR correction (if used) |
| `bandEdgesVec` | ISO band edge vector [Hz] |
| `fTransition` | Schroeder frequency used [Hz] |
| `fs`, `filename`, `bandMode`, `V_room` | Run metadata |

### Figures produced

| Figure | Contents |
|--------|----------|
| *acoustic_indices — aligned RIR* | Original vs onset-aligned impulse response |
| *acoustic_indices — decay and clarity* | Per-band T30 and C80 |
| *polymaxRIR — Spectrum* | Target vs synthetic spectrum, spectral residual, time-domain comparison with nRMSE |
| *polymaxRIR — Per-band statistics* | Mode count per band; mean ± 1σ T60; mean ± 1σ inter-mode spacing |
| *polymaxRIR — LF T60 comparison* | PolyMax mean ± 1σ T60 vs acoustic-indices 2×T30 for all LF bands, annotated with pole counts |

---

## Standalone usage

### `polymax.m`

```matlab
% Multi-order stability diagram, NoCluster = 1 (poles returned directly)
[fBand, cBand, ~] = polymax(H, fv, 2, 200, 5, 0.05, 0.05, 0, 10, 0, 1, 0);
%                           |    |   |   |   |   |     |   |  |   |  |   |
%                       fftinput fv Nmin Nmax step tolf told Nc T60cut AC NC LP
```

Arguments: `fftinput, fv, NorderMin, NorderMax, orderstep, tolf, told, Ncluster, T60cut, AutoCluster, NoCluster, LivePlot`.

### `acoustic_indices.m`

```matlab
[h, fs] = audioread('myRIR.wav');
res = acoustic_indices(h, fs, 3, [20 8000], 0.08, -20);

fprintf('T60 at 1 kHz : %.2f s\n', ...
    interp1(res.bandCentersHz, res.T30, 1000, 'linear', 'extrap'));
fprintf('BR = %.3f  |  TR = %.3f\n', res.BR, res.TR);
```

---

## Algorithm notes

### PolyMax in brief

At each polynomial order N, a scalar right-matrix-fraction model is fitted to the measured FRF:

```
H(z) ≈ B(z) / A(z),    z = exp(−jω·Δt)
```

The denominator coefficients of A(z) are recovered from the Schur complement of the normal equations, stabilised with amplitude-invariant Tikhonov regularisation (the regularisation parameter scales with the mean diagonal of the normal-equation matrices, so results are independent of signal level). Eigenvalues of the companion matrix give the discrete-time poles; only upper-half-plane, left-half-plane poles are retained as physical.

A pole is accepted into the stability diagram when both its relative frequency shift and relative damping shift from one order to the next are below the tolerances `tolf` and `told`.

### Amplitude solve

Given poles (fₚ, cₚ), the modal partial-fraction model is:

```
H(jω) = Σₚ [ rₚ/(jω − λₚ) + rₚ*/(jω − λₚ*) ],    λₚ = −cₚ + j√(ωₚ² − cₚ²)
```

The complex residues rₚ = aₚ + j·bₚ are found by least squares over a frequency window for each band (with a 10 % overlap on each side to reduce edge artefacts). Tikhonov regularisation with `ε = 10⁻⁸ × mean(diag(R'R))` is applied.

### Damping homogenisation

After the PolyMax stability diagram, per-pole damping estimates within each band are replaced by their band mean before the amplitude solve. This prevents a handful of poorly-conditioned poles from producing outlier residues. The original per-pole spread is preserved separately for the LF T60 comparison figure.

### T60 / damping duality

Throughout the codebase, damping is stored as c [rad/s]:

```
T60 [s]  =  3·ln(10) / c  ≈  6.908 / c
```

### ISO band edges

Bands follow ISO 266 / IEC 61260 exactly. Centre frequencies are:

```
fc(n) = 1000 · 2^(n/B)   Hz,    B = bandsPerOctave
```

Band edges are the exact geometric means `fc(n) / 2^(1/2B)` and `fc(n) · 2^(1/2B)`. The user's `fLow`/`fHigh` values select which ISO bands to include — they never clip or adjust band edges.

---

## References

Peeters B., Van der Auweraer H., Guillaume P., Leuridan J. (2004).
*The PolyMAX frequency-domain method: a new standard for modal parameter estimation?*
Shock and Vibration, **11**(3–4), 395–409.

Schroeder M. R. (1965).
*New method of measuring reverberation time.*
Journal of the Acoustical Society of America, **37**(3), 409–412.

ISO 266:1997. *Acoustics — Preferred frequencies.*

IEC 61260-1:2014. *Electroacoustics — Octave-band and fractional-octave-band filters.*
