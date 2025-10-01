function [fitresult, lower_thres, upper_thres] = KLS_fit_loglog1(X, Y, cutoff, plot_flag)
%KLS_FIT_LOGLOG1 Fit a single log–logistic distribution to a PDF histogram.
%
%   [FITRESULT, LOWER_THRES, UPPER_THRES] = KLS_FIT_LOGLOG1(X, Y, CUTOFF, PLOT_FLAG)
%   fits a single log–logistic distribution to histogrammed data provided
%   as bin centers X (ascending) and PDF sample heights Y (same size as X).
%   The function minimizes a discrete negative log-likelihood (NLL) using
%   fmincon (Optimization Toolbox), assuming approximately uniform bin width.
%
%   Inputs
%   ------
%   X         : Bin centers (vector, strictly increasing). Only positive support
%               is modeled; nonpositive X (and matching Y) are dropped with a warning.
%   Y         : Bin heights as a PDF (vector, same size as X). The total mass is
%               validated and (if needed) renormalized using sum(Y)*Δx.
%   CUTOFF    : (optional, scalar in (0,1)) Percentile used for the upper threshold
%               from the component's CDF. Default: 0.95.
%   PLOT_FLAG : (optional, 0/1) If 1, produce a plot overlaying the input PDF,
%               fitted PDF, and the threshold x-line. Default: 0 (no plot).
%
%   Outputs
%   -------
%   FITRESULT : Struct with fitted parameters and helpers:
%       .a, .b                           - Scale and shape of the log–logistic
%       .componentSelected               - Always 1 (kept for template parity)
%       .pdf_fun(z,a,b)                  - Handle: single log–logistic PDF
%       .icdf_fun(p,a,b)                 - Handle: single log–logistic ICDF
%       .mix_pdf_fun(z)                  - Handle: fitted PDF (alias for convenience)
%       .gof.NLL, .gof.AIC, .gof.BIC     - Goodness-of-fit metrics
%       .gof.exitflag                    - fmincon exitflag
%   LOWER_THRES : NaN (left tail not used for right-skewed log–logistic thresholding)
%   UPPER_THRES : Scalar threshold at CUTOFF of the fitted CDF
%
%   Modeling & Optimization
%   -----------------------
%   Log–logistic PDF:
%       ll_pdf(z,a,b) = (b./a).*(z./a).^(b-1) ./ (1 + (z./a).^b).^2,  z>0
%   Discrete NLL minimized:
%       NLL = -sum( (Y.*Δx) .* log( ll_pdf(X,a,b) + eps ) )
%   Parameters and bounds:
%       p = [a, b],  a ∈ [1e-6, Inf),  b ∈ [0.1, Inf)
%   Initialization (deterministic):
%       a = median(posX),  b = 2,  where posX = X(X>0)
%
%   GOF Metrics
%   -----------
%   k = 2 parameters. We report:
%       AIC = 2k + 2*NLL
%       BIC = k*log(nEff) + 2*NLL
%   where nEff is the number of bins with positive mass: nEff = sum((Y.*Δx) > 0).
%
%   Notes
%   -----
%   * Requires Optimization Toolbox (fmincon). If missing, an error is thrown.
%   * Validates strictly increasing X and approximate uniform spacing (Δx uniformity).
%   * Validates/renormalizes Y to integrate to 1 using sum(Y)*Δx.
%   * Deterministic: no randomized initialization.
%
%   Example
%   -------
%   See the runnable usage example at the end of this file (guarded).
%
%   Author: KLS-style template for single log–logistic fit
%   Version: 1.0

    % ------------------------
    % Input validation & setup
    % ------------------------
    narginchk(2,4);

    if nargin < 3 || isempty(cutoff) || ~isscalar(cutoff) || ~isfinite(cutoff) || cutoff <= 0 || cutoff >= 1
        cutoff = 0.95;
    end
    if nargin < 4 || isempty(plot_flag)
        plot_flag = 0;
    end

    validateattributes(X, {'numeric'}, {'vector','real','finite','nonempty'}, mfilename, 'X', 1);
    validateattributes(Y, {'numeric'}, {'vector','real','finite','nonempty','numel',numel(X)}, mfilename, 'Y', 2);

    X = X(:);
    Y = Y(:);

    % Enforce strictly increasing X
    if any(diff(X) <= 0)
        error('X must be strictly increasing.');
    end

    % Drop nonpositive X with warning (log–logistic is only defined for z>0)
    nonpos = X <= 0;
    if any(nonpos)
        warning('Nonpositive X detected and dropped (%d bins). Matching Y elements dropped as well.', sum(nonpos));
        X = X(~nonpos);
        Y = Y(~nonpos);
        if numel(X) < 3
            error('Insufficient positive-support bins after dropping nonpositive X.');
        end
    end

    % Approximate uniform spacing check
    dX  = diff(X);
    dx  = median(dX);
    if dx <= 0
        error('Invalid bin spacing.');
    end
    tol = max(1e-6, 1e-3 * dx); % relative tolerance for uniformity
    if any(abs(dX - dx) > tol)
        error('X must have approximately uniform spacing. Max deviation %.3g exceeds tolerance %.3g.', max(abs(dX - dx)), tol);
    end

    % Normalize Y to integrate to 1 using sum(Y)*dx
    mass = sum(Y) * dx;
    if ~(isfinite(mass) && mass > 0)
        error('Y must define a nonzero PDF mass over X.');
    end
    if abs(mass - 1) > 0.05
        warning('Input Y integrates to %.4f over X; renormalizing to 1.', mass);
    end
    Y = Y / mass; % now sum(Y)*dx == 1 (up to numerical precision)

    % Effective number of bins with positive mass (for BIC)
    nEff = sum((Y * dx) > 0);

    % -------------
    % Toolbox check
    % -------------
    if exist('fmincon', 'file') ~= 2
        error(['KLS_fit_loglog1 requires Optimization Toolbox (fmincon). ', ...
               'Please install the toolbox or ensure it is on the MATLAB path.']);
    end

    % ------------------------
    % Model: helpers & handles
    % ------------------------
    ll_pdf  = @(z,a,b) (b./a) .* (z./a).^(b-1) ./ (1 + (z./a).^b).^2 .* (z > 0);
    ll_icdf = @(p,a,b) a .* (p./(1-p)).^(1./b);

    % Discrete NLL (expectation under histogram mass)
    binMass = Y * dx; % column vector, sums to ~1
    nll_fun = @(p) -sum( binMass .* log( max(realmin, ll_pdf(X, p(1), p(2))) + eps ) );

    % -------------------
    % Initialization/Bnds
    % -------------------
    posX = X; % after dropping nonpositive, all positive
    if isempty(posX)
        error('No positive support available in X for log–logistic fitting.');
    end
    medx = median(posX);

    % Initial guess (deterministic)
    p0 = [medx, 2]; % [a, b]

    % Bounds
    lb = [1e-6, 0.1];
    ub = [Inf,  Inf ];

    % -----------
    % Optimize
    % -----------
    opts = optimoptions('fmincon', ...
        'Algorithm', 'interior-point', ...
        'Display', 'off', ...
        'MaxFunctionEvaluations', 2e4, ...
        'MaxIterations', 2e3, ...
        'OptimalityTolerance', 1e-10, ...
        'StepTolerance', 1e-12, ...
        'ConstraintTolerance', 1e-12);

    % Ensure initial guess within bounds
    p0 = min(max(p0, lb), ub);

    [p_hat, fval, exitflag] = fmincon(nll_fun, p0, [], [], [], [], lb, ub, [], opts);

    a = p_hat(1); 
    b = p_hat(2);

    % ----------------
    % Thresholding
    % ----------------
    upper_thres = ll_icdf(cutoff, a, b);
    lower_thres = ll_icdf(0.05, a, b);

    % ----------------
    % Goodness-of-fit
    % ----------------
    NLL = fval;
    k   = 2; % parameters
    AIC = 2*k + 2*NLL;
    BIC = k*log(max(1, nEff)) + 2*NLL;

    % ----------------
    % Assemble result
    % ----------------
    fitresult = struct();
    fitresult.a = a;
    fitresult.b = b;
    fitresult.componentSelected = 1; % retained for interface parity
    fitresult.pdf_fun  = @(z,varargin) ll_pdf(z, a_if(varargin,1,a), b_if(varargin,2,b));
    fitresult.icdf_fun = @(p,varargin) ll_icdf(p, a_if(varargin,1,a), b_if(varargin,2,b));
    fitresult.mix_pdf_fun = @(z) ll_pdf(z, a, b); % alias for convenience
    fitresult.gof = struct('NLL', NLL, 'AIC', AIC, 'BIC', BIC, 'exitflag', exitflag);

    % -------
    % Plotting
    % -------
    if plot_flag == 1
        zfine = linspace(min(X), max(X), 2000).';
        pdf_fine = ll_pdf(zfine, a, b);

        figure('Color','w');
        hold on;
        plot(X, Y, 'o', 'MarkerSize', 4, 'LineWidth', 1, 'DisplayName', 'Input PDF');
        plot(zfine, pdf_fine, '-', 'LineWidth', 2, 'DisplayName', 'Log–logistic (fit)');

        lbl = sprintf('%.1f%% CDF', cutoff*100);
        xline(upper_thres, '--', lbl, 'LabelVerticalAlignment','top', 'HandleVisibility','on','DisplayName',lbl);
        %lbl = sprintf('%.1f%% CDF', 0.05*100);
        %xline(lower_thres, '--', lbl, 'LabelVerticalAlignment','top', 'HandleVisibility','on','DisplayName',lbl);

        grid on; box on;
        xlabel('X');
        ylabel('Probability density');
        title('Single Log–Logistic Fit');
        legend('Location','best');
        hold off;
    end
end

% --- small local helpers to allow optional (z,a,b) handles above without varargin errors
function v = a_if(vin,idx,default)
    if numel(vin)>=idx && ~isempty(vin{idx})
        v = vin{idx}; 
    else
        v = default; 
    end
end
function v = b_if(vin,idx,default)
    if numel(vin)>=idx && ~isempty(vin{idx})
        v = vin{idx}; 
    else
        v = default; 
    end
end
%{
% -----------------------------
% Runnable usage example (guard)
% -----------------------------
if ~isdeployed && ~isfolder('dummy')
    try
        % Synthetic grid and true parameters
        X = linspace(1e-3, 20, 300).';
        true_a = 3.5; 
        true_b = 2.2;

        ll_pdf = @(z,a,b) (b./a).*(z./a).^(b-1) ./ (1 + (z./a).^b).^2 .* (z > 0);
        Y = ll_pdf(X, true_a, true_b);

        % Normalize to PDF (guard against numerical drift)
        dx = median(diff(X));
        Y = Y / (sum(Y)*dx);

        cutoff = 0.95;
        plot_flag = 1;
        [fitresult, ~, upper_thres] = KLS_fit_loglog1(X, Y, cutoff, plot_flag);

        fprintf('Fitted parameters:\n');
        fprintf('  a  = %.6f, b = %.6f\n', fitresult.a, fitresult.b);
        fprintf('  Upper threshold (%.1f%%): %.6f\n', cutoff*100, upper_thres);
        fprintf('  GOF: NLL=%.6f, AIC=%.6f, BIC=%.6f, exitflag=%d\n', ...
            fitresult.gof.NLL, fitresult.gof.AIC, fitresult.gof.BIC, fitresult.gof.exitflag);
    catch ME
        warning('Example run failed: %s', ME.message);
    end
end
%}