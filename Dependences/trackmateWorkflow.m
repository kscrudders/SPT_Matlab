function results = trackmateWorkflow(img_stack, save_dir, options)
%TRACKMATEWORKFLOW  Minimal TrackMate‑like detection & tracking pipeline.
%
% SYNTAX
%   results = trackmateWorkflow(img_stack, save_dir, 'Param', value, ...)
%
% REQUIRED INPUTS
%   img_stack : H×W×T numeric array already in memory.
%   save_dir  : directory path where .mat/.csv outputs will be written.
%
% OPTIONAL NAME‑VALUE PAIRS
%   'Detector'      : struct with parameters for detectSpots().
%   'Tracker'       : struct with parameters for trackSpots().
%   'ShowPlots'     : logical master switch (default true).
%   'Show'          : struct with logical fields
%                       .qualityHist         (default true)
%                       .detectionOverlay    (default true)
%                       .trackProjection     (default true)
%   'FrameToShow'   : scalar frame index for detection overlay (default ⌊T/2⌋).
%
% OUTPUT
%   results : struct with fields
%               .meta   – image & parameter info
%               .spots  – table of detections (one per row)
%               .tracks – structure array of trajectories
%
% -------------------------------------------------------------------------
% Kevin Scrudders – July 2025
% -------------------------------------------------------------------------
arguments
    img_stack (:,:,:) {mustBeNumeric};
    save_dir  (1,:) {mustBeText}

    options.Detector   struct = struct()
    options.Tracker    struct = struct()
    options.ShowPlots  (1,1) logical = true

    options.Show struct = struct()

    options.FrameToShow (1,1) double = floor(size(img_stack, 3) / 2)
end

if isstring(save_dir)
    save_dir = convertStringsToChars(save_dir);
end

% inputs
Detector = options.Detector;
Tracker = options.Tracker;
    Tracker.pixelsize = Detector.pixelsize; % camera pixel size, care over from detector settings
if ~isempty(Tracker.expected_diffusion) && ~isempty(Tracker.imaging_interval)
    % assume 2D brownian diffusion follows into a Rayleigh distribution of
    % stepsizes, let's pull the expected stepsize for a 95th percentile
    % step

    Tracker.maxDisp = sqrt(2*Tracker.expected_diffusion * Tracker.imaging_interval) * sqrt(-2* log(1-0.95)); % max step in µm/frame
end

default_Show = struct( ...
    'qualityHist',false, ...
    'detectionOverlay',false, ...
    'trackProjection',false, ...
    'filterSpots_SBRHist',false, ...
    'filterSpots_peakHist',false, ...
    'filterSpots_LoGHist',false, ...
    'filterSpots_PCAstats',false);
options.Show = mergeDefaultStructures(options.Show,default_Show);


ShowPlots = options.ShowPlots;
    ShowOpts = options.Show;
FrameToShow = options.FrameToShow;

% outputs
results = struct();
results.FrameToShow = FrameToShow;

% Normalise display flags
if ~ShowPlots
    ShowOpts = structfun(@(~)false, ShowOpts, 'UniformOutput', false);
end

% Select frame to overlay detections
if isnan(results.FrameToShow)
    results.FrameToShow = round(size(img_stack,3)/2);
    disp(['showing frame ' num2str(results.FrameToShow)])
end
if results.FrameToShow > size(img_stack,3)
    results.FrameToShow = size(img_stack,3);
    disp(['showing frame ' num2str(results.FrameToShow)])
end
if results.FrameToShow < 1
    results.FrameToShow = 1;
    disp(['showing frame ' num2str(results.FrameToShow)])
end

% Ensure save_dir exists
if exist(save_dir, 'dir')~=7
    mkdir(save_dir);
end

% ---------------------------------------------------------------------------
% ---- Spot detection ---------------------------------------------------
% ---------------------------------------------------------------------------
[spotsTable, Detector] = detectSpots(img_stack, Detector, ShowPlots, save_dir, options);
%disp('localization done')

% ---------------------------------------------------------------------------
% ---- Spot filter ---------------------------------------------------
% ---------------------------------------------------------------------------
[spotsTable, Detector] = filterSpots(img_stack, spotsTable, Detector, ShowPlots, save_dir, options);
%disp('spot filtering done')

if isempty(spotsTable)
    disp('No Spots Remain after filtering.')
    return;
end
% ---------------------------------------------------------------------------
% ---- Tracking -----------------------------------------------------------
% ---------------------------------------------------------------------------
[tracks, trackerParams, spotsTable] = trackSpots(spotsTable, Tracker);
%disp('tracking done')

% ---------------------------------------------------------------------------
% ---- Track filter ---------------------------------------------------
% ---------------------------------------------------------------------------
[tracks_final, trackerParams_final, spotsTable_final, TrackFilterMetaData] = filterTracks(img_stack, tracks, trackerParams, spotsTable, ShowPlots, save_dir, options, Tracker);
%disp('track filtering done')

% ---------------------------------------------------------------------------
% ---- Figures -----------------------------------------------------------
% ---------------------------------------------------------------------------
if ShowPlots && ShowOpts.trackProjection
    plotTracksProjection(img_stack, tracks_final);
    savefig(fullfile(save_dir, 'MaxProjection_TrackOverlay'))
end

if ShowPlots && ShowOpts.trackVideo
    %p.addParameter('MarkerSize', 5);
    %p.addParameter('LineWidth', 1.5);
    %p.addParameter('Colors', []);
    %p.addParameter('Codec', 'MPEG-4');
    %p.addParameter('PixelSize', 0.157);
    %p.addParameter('ScaleBar', 10);
    %p.addParameter('Figure', []);
    STLN_tracks = KLS_tracks_2_STLN(img_stack, tracks_final);

    videoTracks(img_stack, KLS_Fill_TrackZeros(STLN_tracks), fullfile(save_dir, 'tracks.mp4'), ...
        'FPS', 15, 'Tail', 50, 'Fade', 1);
end

%% ---- Assemble results --------------------------------------------------

results.meta       = struct( ...
        'imgSize',       size(img_stack), ...
        'Detector',      Detector, ...
        'Tracker',       trackerParams_final);
results.spots      = spotsTable_final;
results.tracks     = tracks_final;
results.NT2_tracks = KLS_tracks_2_STLN(img_stack, tracks_final);
results.TrackFilterMetaData = TrackFilterMetaData;

results.pre_filter.spots = spotsTable;
results.pre_filter.tracks = tracks;

%%---- Save outputs ------------------------------------------------------
matFile = fullfile(save_dir,'results.mat');
save(matFile, 'results', '-v7.3');

spotsCSV  = fullfile(save_dir,'spots.csv');
writetable(spotsTable_final, spotsCSV);

if ShowPlots
    fprintf('[tracking complete] Results saved to %s\n', save_dir);
end

end


function S = mergeDefaultStructures(S,D)
    if isempty(S), S = struct(); end
    f = fieldnames(D);
    for k = 1:numel(f)
        if ~isfield(S,f{k}) || isempty(S.(f{k}))
            S.(f{k}) = D.(f{k});
        end
    end
end