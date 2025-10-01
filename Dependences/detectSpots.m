function [spots, Detector] = detectSpots(img_stack, Detector, plotFlag, save_dir, options)
%DETECTSPOTS  LoG‑based spot detection with optional MLE refinement.
%
%   [spots, DetectorParams] = detectSpots(img_stack, DetectorStruct, plotFlag)
%
% INPUTS
%   img_stack  : H×W×T numeric array held in memory.
%   Detector   : struct with at least the following fields:
%       diameter          – expected object diameter (pixels, scalar, req.)
%       qualityThreshold  – scalar threshold or [] to auto‑estimate (default [])
%       minIntensity      – lower raw‑intensity cutoff (optional, [])
%       maxIntensity      – upper raw‑intensity cutoff (optional, [])
%       subpixel          – true | false (default true)  perform MLE refinement
%   plotFlag   : true  => show histogram + overlay figure (default = false)
%
% OUTPUTS
%   spots          – MATLAB table, one row per detection (SpotID, Frame, x_px, y_px, Quality …)
%   DetectorParams – struct capturing *all* numeric settings actually used.  Store this
%                    alongside results for future reproducibility / debugging.
%
% NOTES
%   • Requires Image Processing + Optimization toolboxes when subpixel refinement is enabled.
%   • Quality metric = −LoG response at the candidate location (larger ⇒ better).
%
% Kevin Scrudders, 2025‑07‑13

% -------------------------------------------------------------------------
% defaults & sanity checks
if nargin < 3 || isempty(plotFlag)
    plotFlag = false; 
end

assert(ndims(img_stack) == 3, 'img_stack must be H×W×T');

% fill in missing Detector fields
if ~isfield(Detector, 'qualityThreshold')
    Detector.qualityThreshold = []; % empty sets up auto thresholding 
end
if ~isfield(Detector, 'minIntensity')
    Detector.minIntensity = 0;  
end
if ~isfield(Detector, 'maxIntensity')
    Detector.maxIntensity = 2^16;  
end
if ~isfield(Detector, 'subpixel')
    Detector.subpixel = true;    
end
if ~isfield(Detector, 'diameter')
    Detector.diameter = 3; % 3 pixel diameter default    
end

[H, W, T] = size(img_stack);

% -------------------------------------------------------------------------
% 1) build LoG filter
sigma      = (Detector.diameter ./ Detector.pixelsize) / 2*sqrt(2*log(2)); % diameter = FWHM in pixel, convert to sigma, assuming gaussian shape
filterSz   = (3*ceil(sigma))+1;                        % ensure odd size
logKernel  = fspecial('log', filterSz, sigma);

% pre‑allocate containers
maxSpotsPerFrame = ceil( (H*W) / ((Detector.diameter ./ Detector.pixelsize)^2) );
spotFrame  = zeros(maxSpotsPerFrame*T,1,'uint32');
spotXY     = zeros(maxSpotsPerFrame*T,2);
spotQual   = zeros(maxSpotsPerFrame*T,1);
count = 0;

% -------------------------------------------------------------------------
% 2a) apply the filter
L= imfilter(img_stack, logKernel, 'replicate', 'conv'); % LoG response
Q = -L; % quality metric (peaks are negative)

% 2-D local maxima on every slice at once
se    = ones(3,3,1);
isMax = Q == imdilate(Q, se);                 % candidates
%isMax = isMax & (Q > imerode(Q, se));         % kill plateaus if desired

% Gather results
[idxR, idxC, idxT] = ind2sub(size(isMax), find(isMax));
spotFrame = idxT;
spotXY    = [idxC idxR];  % x=col, y=row
spotQual  = Q(isMax);

% -------------------------------------------------------------------------
% 3) quality thresholding (auto if needed)
if isempty(Detector.qualityThreshold)
    % robust estimate: fit normal to central bulk and set thr = mu + 3*sd
    spotQual_fit = spotQual(spotQual >0); % filter 0 spot quality
    [Y, binEdges] = histcounts(spotQual_fit, 'BinMethod', 'fd','Normalization','PDF'); % The Freedman-Diaconis rule is less sensitive to outliers in the data, and might be more suitable for data with heavy-tailed distributions. It uses a bin width of 2*iqr(X(:))*numel(X)^(-1/3).

    X = movmean(binEdges,2);
    X = X(2:end);

    cutoff     = [0.98]; 
    % Log-logistic distribution
    [~, ~, qualityThr] = KLS_fit_loglog1(X, Y, cutoff, 0);

    autoThr = true;
else
    qualityThr = Detector.qualityThreshold;
    autoThr    = false;
end
%
% -------------------------------------------------------------------------
% [Plot] Quality spot dection threshold
if plotFlag && options.Show.qualityHist == 1
    [~, ~, ~] = KLS_fit_loglog1(X, Y, cutoff, plotFlag);
    
    % Label axes and add title
    xlabel('LoG Quality');
    ylabel('PDF');
    
    axis tight; % Remove extra whitespace

    savefig(fullfile(save_dir, 'LoG_Quality_Threshold'))
end
% -------------------------------------------------------------------------
% remove spots of low quality

keep = spotQual >= qualityThr;
spotFrame = spotFrame(keep);
spotXY    = spotXY(keep,:);
spotQual  = spotQual(keep);

% -------------------------------------------------------------------------
% 4) MLE sub‑pixel refinement
win = 7;
if Detector.subpixel && ~isempty(spotQual)
    [spotXY, spotAmplitude, fitInfo, keep_all_positive] = refineGaussianGPU(img_stack, spotXY, spotFrame, sigma, win);
end
spotFrame = spotFrame(keep_all_positive);
spotQual  = spotQual(keep_all_positive);

% -------------------------------------------------------------------------
% 5) build output table
N = numel(spotQual);
spots = table( ...
    (1:N)', ...
    spotFrame, ...
    spotXY(:,1), ...
    spotXY(:,2), ...
    spotQual, ...
    fitInfo.amp, ...                   % peak amplitude (α)
    fitInfo.sigma, ...                 % fitted sigma
    fitInfo.bg, ...                     % background per pixel (β)
    fitInfo.spot_NetSignal, ...         % net signal (continuous)
    fitInfo.spotSignal, ...             % integrated continuous signal
    fitInfo.spotBkgd, ...               % integrated continuous background
    fitInfo.integrated_signal_roi, ...  % raw ROI-sum integrated signal
    'VariableNames', { ...
        'SpotID', ...
        'Frame', ...
        'x_px', ...
        'y_px', ...
        'Quality_LoG', ...
        'Peak', ...
        'Fit_sigma', ...
        'Background_px', ...
        'Net_Signal', ...
        'Integrate_fit_Signal', ...
        'Integrate_fit_Bkgd', ...
        'Integrate_Signal_ROI' ...
    } ...
);

% -------------------------------------------------------------------------
% 6) [Plot] Spot Overlay
if plotFlag && options.Show.detectionOverlay == 1
    figure()
    img = img_stack(:,:,options.FrameToShow);
    imshow(img,[]);
    
    hold on
        row_idx = spots{:,'Frame'} == options.FrameToShow;
        % Pull out coordinates and orient as N by T
        x = spots{row_idx,'x_px'};   % N by T
        y = spots{row_idx,'y_px'};   % N by T
        
        % Drop tracks that are entirely NaN (never detected)
        valid = any(~isnan(x),2);            % Nx1 logical
        x = x(valid,:);
        y = y(valid,:);
        
        x = x.'; % reshape the matix so its T by N so plot reads the points correctly
        y = y.';
        
        % plot() accepts a matrix: each column (now becuase we reshaped the vector this is N) becomes one line
        scatter(x, y, 'o', 'LineWidth',2, 'Color', [0 0 1 0.4]); % one color
        %plot(x, y, 'LineWidth',2); % different colors per track
    hold off
    
    axis tight; % Remove extra whitespace

    title('Pre-Localization Filter')

    savefig(fullfile(save_dir, 'Localization Check - Pre-Localization Filter'))
end

% -------------------------------------------------------------------------
% 7) assemble DetectorParams for reproducibility

% store kernel because sigma alone is not enough when debugging filtering artefacts
% (downsides: memory, but kernel is small)

% run‑time / data info
s = warning('off','all'); %#ok suppress during struct display
    Detector.dateTime      = datestr(now, 'yyyy‑mm‑dd HH:MM:SS');
    Detector.imageSize     = [H W T];
    Detector.nSpots        = N;
warning(s);

% settings
paramFields = {'diameter','qualityThreshold','minIntensity','maxIntensity','subpixel'};
for f = paramFields
    if isfield(Detector, f{1})
        Detector.(f{1}) = Detector.(f{1});
    else
        Detector.(f{1}) = [];
    end
end

Detector.sigma           = sigma;
Detector.filterKernel    = logKernel;
Detector.filterSize      = filterSz;
Detector.autoThreshold   = autoThr;
Detector.qualityThrUsed  = qualityThr;

end


%%


