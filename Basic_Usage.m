data_dir = 'C:\your data's folder';

img_stack = data; % You'll have to import the data into matlab (see KLS_ND2ImportAll in my other repo)

% Output directory
save_dir = [data_dir, 'Matlab_Tracking'];

if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

% Detector settings
detectorParams = struct( ...
    'diameter', 0.457, ... % microns, FWHM of PSF if you have it, otherwise guess - 
    ... % double lateral resolution - https://www.microscope.healthcare.nikon.com/microtools/resolution-calculator/
    'qualityThreshold', [], ...    % auto-estimate
    'minIntensity', 100, ...
    'maxIntensity', 30000, ...
    'subpixel', true, ... % fit spots with a 2D guassain
    'pixelsize', 0.157, ... % microns/pixel
    'minSNR', 1, ... % minimum spot SNR 
    'maxSNR',30, ... % maximum spot SNR
    'maxSigma', 3, ... % maximum spot sigma (px)
    'minSigma', 1 ... % minimum spot sigma (px)
    );

% Tracker settings
trackerParams = struct( ...
    'maxDisp',  1.5, ... microns/frame
    'maxGap',   20, ... % frames
    'expected_diffusion', 1.0, ... % micron/s^2
    'imaging_interval', 0.050, ... % frame to frame time (s)
    'minLen', 25 ... % minimum number of spots for a valid tracke
    );

% Display options
showPlots = true;
showOpts = struct( ...
    'qualityHist',      false, ...
    'detectionOverlay', false, ...
    'trackProjection',  false, ...
    'trackVideo', false);

% Frame to overlay detection results
frameToShow = 1;

% Run the workflow

results = trackmateWorkflow(img_stack, save_dir, ...
    'Detector',    detectorParams, ...
    'Tracker',     trackerParams, ...
    'ShowPlots',   showPlots, ...
    'Show',        showOpts, ...
    'FrameToShow', frameToShow);
