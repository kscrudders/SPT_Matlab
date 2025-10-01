function [tracks, Tracker, spotsTblOut, TrackFilterMetaData] = filterTracks(img_stack, tracks, trackerParams , spotsTableIn, ShowPlots, save_dir, options, Tracker)
%FILTERTRACKSBYLENGTH  Remove short tracks and keep table/array in sync.
%
%   [tracks, spotsTblOut, XYout, keepMask, info] = ...
%       filterTracksByLength(tracks, spotsTableIn, minLen)
%
% INPUTS
%   tracks   : struct array from trackSpots() with fields
%                .id, .frames, .x, .y, ...  (ID values are positive integers)
%   spotsTableIn : table that includes a numeric 'TrackID' column
%   STLN_tracks         : N×T×2 numeric array (x,y) with NaNs for gaps (N = #tracks)
%   minLen     : minimum number of valid (x,y) samples required to keep a track
%
% OUTPUTS
%   tracks  : filtered track struct array (IDs reindexed to 1..N_kept)
%   spotsTblOut: filtered spots table, with TrackID remapped to new IDs
%   XYout      : filtered N_kept×T×2 array (rows correspond to tracks.id)

%   info       : struct with bookkeeping (removedIDs, idMapOld2New, lengths, etc.)
%
% NOTES
%   • A frame contributes to length only if BOTH x and y are finite.
%   • Tracks are reindexed to be contiguous after filtering. Use info.idMapOld2New
%     if you need to translate original IDs to new IDs.
%   • Robust to cases where STLN_tracks has extra empty rows; falls back to struct lengths
%     if needed.
%
% Kevin Scrudders — Aug 2025

    % ---------- Validate inputs ----------
    arguments
        img_stack (:,:,:) {mustBeNumeric}
        tracks   (1,:) struct
        trackerParams struct = struct()
        spotsTableIn table = table()
        ShowPlots logical = false;
        save_dir (1,:) = {mustBeText}
        options struct = struct()
        Tracker struct = struct()
    end

    if ~ismember('TrackID', spotsTableIn.Properties.VariableNames)
        error('spotsTableIn must contain a numeric ''TrackID'' column.');
    end
    if ~isnumeric(spotsTableIn.TrackID)
        error('spotsTableIn.TrackID must be numeric.');
    end

    if ~isfield(Tracker,'minLen')
        Tracker.minLen = 5;
    else
        if isempty(Tracker.minLen) | Tracker.minLen < 0
            Tracker.minLen = 1;
        end
    end

    minLen = Tracker.minLen;

    % ---------- Basic sizes / IDs ----------
    STLN_tracks = KLS_tracks_2_STLN(img_stack, tracks);
    Nxy = size(STLN_tracks,1);
    T   = size(STLN_tracks,2); %#ok<NASGU>  % (may be handy downstream)
    ids = reshape([tracks.id], 1, []);           % original track IDs
    if any(ids < 1) || any(~isfinite(ids))
        error('tracks.id must be positive finite integers.');
    end
    maxID = max(ids);

    % ---------- Compute per-track lengths ----------
    % Preferred: count frames where BOTH x,y are finite in STLN_tracks
    lenXY_all = sum(all(isfinite(STLN_tracks),3), 2);   % N×1, length of tracks
    ids = (1:length(lenXY_all))';

    % ---------- Decide which tracks to keep ----------
    keepMaskTracks = (lenXY_all >= minLen);   % over the ORDER of tracks
    keptIDs        = ids(keepMaskTracks);
    removedIDs     = ids(~keepMaskTracks);

    % ---------- Build old→new ID map for remapping table & struct ----------
    idMapOld2New = zeros(max(maxID, Nxy), 1, 'uint32');  % 0 = removed
    idMapOld2New(keptIDs) = uint32(1:numel(keptIDs)); %

    % ---------- Filter / reindex the track struct ----------
    tracks = tracks(keepMaskTracks);
    new_ids = idMapOld2New(keepMaskTracks);
    C = num2cell(new_ids);
    [tracks.id] = C{:};

    % ---------- Filter & remap spotsTbl ----------
    % 1) Drop rows from removed tracks
    inKept = ismember(spotsTableIn.TrackID, keptIDs);
    spotsTblOut = spotsTableIn(inKept, :);

    % 2) Remap TrackID to the new contiguous IDs
    oldIDs_vec = spotsTblOut.TrackID;
    % Use safe mapping via array indexing (removedIDs are already dropped)
    newIDs_vec = idMapOld2New(oldIDs_vec);
    spotsTblOut.TrackID = double(newIDs_vec);

    % ---------- Info / bookkeeping ----------
    TrackFilterMetaData = struct();
    TrackFilterMetaData.minLen          = minLen;
    TrackFilterMetaData.lengthsByTrack  = table(ids(:), keepMaskTracks(:), ...
                                 'VariableNames', {'OldID','KeepMask'});
    TrackFilterMetaData.removedIDs      = removedIDs(:);
    TrackFilterMetaData.keptIDs         = keptIDs(:);
    TrackFilterMetaData.idMapOld2New    = idMapOld2New;
    TrackFilterMetaData.N_in            = numel(ids);
    TrackFilterMetaData.N_out           = numel(keptIDs);
end
