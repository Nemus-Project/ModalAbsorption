function [fBand, cBand, Ncluster] = polymax(fftinput, fv, NorderMin, NorderMax, ...
    orderstep, tolf, told, Ncluster, T60cut, AutoCluster, NoCluster, LivePlot)
%POLYMAX  Frequency-domain PolyMax modal parameter estimator.
%
%   [FBAND, CBAND, NCLUSTER] = POLYMAX(FFTINPUT, FV, NORDERMIN, NORDERMAX,
%       ORDERSTEP, TOLF, TOLD, NCLUSTER, T60CUT, AUTOCLUSTER, NOCLUSTER, LIVEPLOT)
%
%   Estimates modal frequencies and damping coefficients from a single-sided
%   complex frequency response FFTINPUT, evaluated at the frequency vector FV
%   (Hz).  The algorithm fits a right matrix-fraction polynomial model over a
%   range of polynomial orders [NORDERMIN : ORDERSTEP : NORDERMAX] and
%   identifies physical poles through a stability diagram: a pole is declared
%   stable when both its relative frequency shift and relative damping shift
%   between consecutive orders fall within the tolerances TOLF and TOLD.
%   Only poles with T60 = 3*ln(10)/c < T60CUT are retained.
%
%   Pole clustering is controlled by AUTOCLUSTER and NOCLUSTER:
%     NoCluster = 1  :  stable poles from the stability diagram are returned
%                       directly (no clustering).  Used for LF bands where
%                       individual modes are resolved.
%     NoCluster = 0  :  all stable poles collected across all orders are
%                       hierarchically clustered into NCLUSTER representatives.
%                       Used for HF bands where modal overlap prevents direct
%                       resolution.
%     AutoCluster = 1:  NCLUSTER is overridden by round(mean(Npoles(2:end))).
%
%   INPUTS
%     fftinput   : complex single-sided FFT of the RIR [Nfft x 1]
%     fv         : frequency vector corresponding to FFTINPUT [Hz]
%     NorderMin  : minimum polynomial order (>= 2)
%     NorderMax  : maximum polynomial order
%     orderstep  : order increment
%     tolf       : relative frequency tolerance for stability (e.g. 0.01)
%     told       : relative damping  tolerance for stability (e.g. 0.01)
%     Ncluster   : number of clusters (used when AutoCluster=0, NoCluster=0)
%     T60cut     : reject poles with T60 > T60cut [s]
%     AutoCluster: 1 = override Ncluster automatically, 0 = use supplied value
%     NoCluster  : 1 = skip clustering and return poles directly
%     LivePlot   : 1 = draw the stability diagram in the current figure
%
%   OUTPUTS
%     fBand    : modal frequencies [Hz], sorted ascending
%     cBand    : modal damping coefficients [rad/s]  (T60 = 3*ln10/c)
%     Ncluster : number of clusters actually used
%
%   DEPENDENCIES
%     Statistics and Machine Learning Toolbox  (linkage, cluster) when
%     NoCluster = 0.
%
%   NOTES
%     - The frequency vector is internally shifted to start at zero to
%       improve numerical conditioning of the Vandermonde-like matrices.
%     - Tikhonov regularisation is amplitude-invariant: the regularisation
%       parameter scales with the mean diagonal of the normal-equation
%       matrices, so results are independent of input signal level.
%     - Only upper-half-plane poles (Im(lambda) > 0) with negative real
%       part (physical damping) are retained.
%
%   REFERENCES
%     Peeters, B., Van der Auweraer, H., Guillaume, P., & Leuridan, J. (2004).
%     The PolyMAX frequency-domain method: a new standard for modal parameter
%     estimation?  Shock and Vibration, 11(3-4), 395-409.
%
%   See also: polymaxRIR, acoustic_indices

% =========================================================================
%  Initialisation
% =========================================================================

%-- Frequency-domain setup
omv    = 2*pi*fv(:);
fstart = fv(1);
fend   = fv(end);
omv    = omv - omv(1);          % shift to zero for numerical conditioning
Lfv    = numel(fv);

%-- Magnitude spectrum (for plotting and LineLevels)
spc      = abs(fftinput(:));
spc_dB   = 20*log10(spc + eps);
maxSpc   = max(spc_dB);
minSpc   = min(spc_dB);

%-- Stability-diagram book-keeping
maxIters   = floor((NorderMax - NorderMin) / orderstep) + 1;
LineLevels = linspace(minSpc, maxSpc, maxIters);
plotInd    = 0;

fTot     = zeros(1e5, 1);
cTot     = zeros(1e5, 1);
cumPoles = 0;
Npoles   = zeros(maxIters, 1);
indPoles = 1;

iter  = 0;
fprev = [];
cprev = [];
fst   = [];
cst   = [];

% =========================================================================
%  Main PolyMax loop — sweep polynomial order
% =========================================================================

for Nopt0 = NorderMin : orderstep : NorderMax

    % Rescale time step as per PolyMax normalisation (Guillaume et al.)
    k = 1 / (fend - fstart) / 2;

    % Build regression matrices X (denominator basis) and Y (numerator coupling)
    X = zeros(Lfv, Nopt0 + 1);
    Y = zeros(Lfv, Nopt0);
    d = zeros(Lfv, 1);
    for n = 1 : Lfv
        z      = exp(-1j * omv(n) * k);   % z-domain sample
        mH     = -fftinput(n);             % negated FRF sample
        for m = 0 : Nopt0
            X(n, m+1) = z^m;
        end
        for m = 0 : Nopt0 - 1
            Y(n, m+1) = mH * z^m;
        end
        d(n) = z^Nopt0 * mH;
    end

    % Normal equations
    rXX = real(X'*X);
    rXY = real(X'*Y);
    rXd = real(X'*d);
    rYY = real(Y'*Y);
    rYX = real(Y'*X);
    rYd = real(Y'*d);

    % Amplitude-invariant Tikhonov regularisation:
    %   eps_r * (mean diagonal) * I  scales with signal power so the solve
    %   is independent of input amplitude.
    eps_r = 1e-6;
    regX  = eps_r * (trace(rXX) / (Nopt0 + 1))     * eye(Nopt0 + 1);
    regY  = eps_r * (trace(rYY) / max(Nopt0, 1))   * eye(Nopt0);

    % Schur-complement solve for denominator coefficients
    al   = (rYY + regY - rYX * ((rXX + regX) \ rXY)) \ ...
           (rYX * ((rXX + regX) \ rXd) - rYd);

    % Build companion matrix and extract eigenvalues (discrete-time poles)
    avec = -al(end - Nopt0 + 1 : end);
    Gm   = [zeros(Nopt0-1, 1), eye(Nopt0-1); avec.'];
    lmd  = eig(Gm);

    % Map to continuous-time and keep physically meaningful poles:
    %   upper half-plane (Im > 0)  and  left half-plane (Re < 0, damped)
    lmd              = -log(lmd) / k;
    lmd(imag(lmd) <= 0) = [];
    lmd(real(lmd) >= 0) = [];

    % Recover frequency and damping, sort by frequency
    cOpt             = -real(lmd);
    omOpt            = sqrt(imag(lmd).^2 + cOpt.^2) + 2*pi*fstart;
    [omOpt, indSort] = sort(omOpt);
    fOpt             = omOpt / (2*pi);
    cOpt             = cOpt(indSort);
    NOpt             = numel(lmd);

    % -----------------------------------------------------------------
    %  Stability selection
    % -----------------------------------------------------------------
    if maxIters == 1
        % Single-order run: return all poles that pass the T60 cut
        fst   = fOpt;
        cst   = cOpt;
        TF    = find(3*log(10) ./ max(cst, eps) > T60cut);
        cst(TF) = [];
        fst(TF) = [];
        indst = numel(fst);

        Npoles(indPoles) = indst;
        fTot = fst;
        cTot = cst;
        cumPoles = indst;

    else
        % Multi-order run: require pole stability across consecutive orders
        indst = 0;

        if iter == 0
            % First pass — seed previous-order pool
            fprev = fOpt;
            cprev = cOpt;
            fst   = fOpt;
            cst   = cOpt;
            indst = numel(cst);

        else
            % Compare current poles against previous-order pool
            if ~isempty(fprev)
                for n = 1 : NOpt
                    % Guard against near-zero denominator
                    fSafe = max(fOpt(n), eps);
                    cSafe = max(cOpt(n), eps);

                    [relFreqShift, idxNearest] = ...
                        min(abs(fOpt(n) - fprev) / fSafe);
                    relDampShift = ...
                        abs(cOpt(n) - cprev(idxNearest)) / cSafe;

                    % Accept pole if stable and within T60cut
                    if relFreqShift <= tolf && relDampShift <= told && ...
                            3*log(10) / max(cprev(idxNearest), eps) < T60cut
                        indst      = indst + 1;
                        fst(indst) = fprev(idxNearest);
                        cst(indst) = cprev(idxNearest);
                    end
                end
            end
            fprev = fOpt;
            cprev = cOpt;
        end

        fst = fst(1:indst);
        cst = cst(1:indst);

        % Secondary T60 filter (belt-and-suspenders)
        TF        = find(3*log(10) ./ max(cst, eps) > T60cut);
        cst(TF)   = [];
        fst(TF)   = [];
        indst     = numel(fst);

        Npoles(indPoles)                         = indst;
        fTot(cumPoles + 1 : cumPoles + indst)   = fst;
        cTot(cumPoles + 1 : cumPoles + indst)   = cst;
        indPoles = indPoles + 1;
        iter     = iter + 1;
        cumPoles = cumPoles + indst;
    end

    % -----------------------------------------------------------------
    %  Live stability diagram (optional)
    % -----------------------------------------------------------------
    if LivePlot == 1
        cla; hold on;
        fill([fv(:); flipud(fv(:))], ...
             [spc_dB(:); (minSpc - 15)*ones(Lfv, 1)], ...
             [0.88 0.88 0.88], 'EdgeColor', 'none', 'FaceAlpha', 0.8);
        plot(fv, spc_dB, 'Color', [0.3 0.3 0.3], 'LineWidth', 0.7);

        if indst > 0
            plotInd = plotInd + 1;
            frac = (Nopt0 - NorderMin) / max(NorderMax - NorderMin, 1);
            col  = [frac, 0.2, 1 - frac];   % blue → red gradient
            plot(fst, LineLevels(min(plotInd, end)) * ones(indst, 1), ...
                'LineStyle', 'none', 'Marker', 'x', ...
                'Color', col, 'MarkerSize', 5, 'LineWidth', 1.2);
        end

        xlim([fv(1), fv(end)]);
        ylim([minSpc - 5, maxSpc + 5]);
        xlabel('$f$ (Hz)', 'Interpreter', 'latex', 'FontSize', 11);
        ylabel('$|\hat{H}|$ (dB)', 'Interpreter', 'latex', 'FontSize', 11);
        title(sprintf('Order %d/%d  |  stable this step: %d  |  total: %d', ...
            Nopt0, NorderMax, indst, cumPoles), 'FontSize', 9);
        drawnow;
    end

end   % order sweep

fTot = fTot(1:cumPoles);
cTot = cTot(1:cumPoles);

% =========================================================================
%  Pole clustering
% =========================================================================

if NoCluster == 1
    % Return stability-diagram poles directly (LF regime)
    [fBand, indSort] = sort(fst);
    cBand            = cst(indSort);
    Ncluster         = numel(cBand);
    return;
end

% Hierarchical clustering of the full fTot pool (HF regime)
if AutoCluster == 1
    if numel(Npoles) > 1
        Ncluster = round(mean(Npoles(2:end)));
    else
        Ncluster = Npoles(1);
    end
end

if numel(fst) < Ncluster
    Ncluster = numel(fst);
end

if Ncluster == 0 || isempty(fTot)
    fBand = [];
    cBand = [];
    return;
end

Z     = linkage(fTot);
clust = cluster(Z, 'Maxclust', Ncluster);

cCentroid = zeros(Ncluster, 1);
fCentroid = zeros(Ncluster, 1);

for m = 1 : Ncluster
    idx_m = (clust == m);
    if any(idx_m)
        cCentroid(m) = mean(cTot(idx_m));
        fCentroid(m) = mean(fTot(idx_m));
    end
end

[fBand, indSort] = sort(fCentroid);
cBand            = cCentroid(indSort);

% =========================================================================
%  Final cluster overlay on live plot
% =========================================================================

if LivePlot == 1
    for n = 1 : Ncluster
        fn = fBand(n);
        pl = line([fn, fn], [minSpc, maxSpc], ...
            'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.2);
        if n > 1
            set(get(get(pl, 'Annotation'), 'LegendInformation'), ...
                'IconDisplayStyle', 'off');
        end
    end
    xlim([fv(1), fv(end)]);
    ylim([minSpc - 10, maxSpc + 10]);
    xlabel('$f$ (Hz)',        'Interpreter', 'latex', 'FontSize', 12);
    ylabel('$|\hat{y}|$ (dB)', 'Interpreter', 'latex', 'FontSize', 12);
    title(sprintf('Stable clusters: %d  |  [%.1f – %.1f Hz]', ...
        Ncluster, fv(1), fv(end)), 'FontSize', 10);
    drawnow;
    pause(2);
    hold off;
end

end   % polymax