function [tracks, trackerParams, spotsTbl] = trackSpots(spotsTbl, Tracker)
%TRACKSPOTS  Simple nearest‑neighbour & gap‑closing tracker for 2‑D spots.
%
%   [tracks, trackerParams] = trackSpots(spotsTbl, trackerStruct)
%
%   INPUTS
%   ------
%   spotsTbl : table as returned by detectSpots() with columns
%                x_px, y_px, Frame, Quality_LoG, SpotID (any extra ignored)
%   tracker  : struct with optional fields (defaults shown):
%                .maxDisp     5  % max displacement (pixels) between frames
%                .maxGap      2  % permissible frame gaps for linking
%                .useGPU      false (reserved / unused in this minimal impl.)
%
%   OUTPUTS
%   -------
%   tracks        : structure array, one element per trajectorys
%     .id         : track integer ID
%     .frames     : [1×N] vector of frame numbers
%     .x, .y      : [1×N] vectors of centroid positions (pixels)
%     .quality    : [1×N] vector of spot quality metrics
%
%   trackerParams : struct of the effective tracking settings (for metadata)
%
%   Notes
%   -----
%   • A greedy nearest‑neighbour LAP (Jonker–Volgenant) assignment is applied per
%     frame pair (t vs t+1) via `matchpairs`. Only links within `maxDisp`
%     are allowed.  Unlinked spots start new trajectories.
%   • Simple gap‑closing: tracks may remain alive for up to `maxGap` empty
%     frames before they are terminated. No branching/merging is handled.
%   • For additional options work consider TrackMate (2022) (https://www.nature.com/articles/s41592-022-01507-1),  
%     of u‑track - Jaqaman et al. (2008).
%   • The Jonker-Volgenant algorithm comes from Yi Cao and some speed
%   improves by Eric Trautmann (https://matlab.mathworks.com/open/fileexchange/v1?id=30838)
%     
% -------------------------------------------------------------------------
% Kevin Scrudders  |  July 2025
% -------------------------------------------------------------------------
    arguments
        spotsTbl table
        Tracker  struct = struct();
    end
    
    % -------------------------- Parameters --------------------------
    % defaults
    trackerParams = struct( ...
        'maxDisp', Tracker.maxDisp, ...
        'maxGap',  Tracker.maxGap, ...
        'timestamp', string(datetime('now','Format','yyyy-MM-dd HH:mm:ss')));

    fns = fieldnames(Tracker);
    for k = 1:numel(fns), trackerParams.(fns{k}) = Tracker.(fns{k}); end
    trackerParams.timestamp = string(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));
    
    maxDisp = (trackerParams.maxDisp ./ trackerParams.pixelsize); % link distance (px)
    th2 = maxDisp^2; % threshold for when we are checking square distances
    maxGap  = trackerParams.maxGap;
    
    % -------------------------- Pre-process detections --------------------------
    spotsTbl = sortrows(spotsTbl, {'Frame'});  % chronological order (one-time)
    
    % Extract columns once (no table reads in loop)
    x  = spotsTbl.x_px;
    y  = spotsTbl.y_px;
    fr = spotsTbl.Frame;
    q  = spotsTbl.Quality_LoG;
    P  = spotsTbl.Peak;
    NS = spotsTbl.Net_Signal;
    IS = spotsTbl.Integrate_fit_Signal;
    IB = spotsTbl.Integrate_fit_Bkgd;
    
    % Frame blocks (avoid find(fr==t) in loop)
    edges = [1; find(diff(fr)~=0)+1; numel(fr)+1];
    AllF  = fr(edges(1:end-1));
    
    nSpots  = numel(x);
    trackID = zeros(nSpots,1,'uint32');
    
    % -------------------------- Active state --------------------------
    activeID    = uint32([]);   % track IDs for active tracks
    activeXY    = zeros(0,2);   % last known (x,y)
    activeFrame = zeros(0,1);   % last frame index
    activeGap   = zeros(0,1);   % gap counter
    nextID      = uint32(1);
    
    % -------------------------- Iterate over frames --------------------------
    for tIdx = 1:numel(AllF)
        t = AllF(tIdx);
        idxThis = edges(tIdx):edges(tIdx+1)-1;
        xyzThis = [x(idxThis), y(idxThis)];
        nActive = numel(activeID);
        nNow    = size(xyzThis,1);
    
        % Prepare containers for this frame's linking
        pairs = zeros(0,2);  % [activeRow, spotRow] in local (active/idxThis) indices
    
        if nActive>0 && nNow>0
            % ---- 1) Reciprocal NN pre-matches (cheap & safe) ----
            [idxT2S, dT2S] = knnsearch(xyzThis,  activeXY, 'K',1);  % track -> spot
            [idxS2T, dS2T] = knnsearch(activeXY, xyzThis,  'K',1);  % spot  -> track

            mask_rnn = ((1:nActive)' == idxS2T(idxT2S)) & (dT2S <= maxDisp);
            pairs_rnn = [find(mask_rnn), idxT2S(mask_rnn)];   % local indices
    
            % Remaining candidates after removing RNN matches
            keepTrack = true(nActive,1); 
            if isempty(pairs_rnn)
                keepTrack(:) = false;
                keepSpot  = true(nNow,1);    
                keepSpot(:)  = false;            
            else
                keepTrack(pairs_rnn(:,1)) = false;
                keepSpot  = true(nNow,1);    
                keepSpot(pairs_rnn(:,2))  = false;                
            end

    
            if any(keepTrack) && any(keepSpot)
                % ---- 2) Candidate neighbourhoods via KD-tree ----
                Mdl    = createns(xyzThis(keepSpot,:), 'NSMethod','kdtree');
                nbrIdx = rangesearch(Mdl, activeXY(keepTrack,:), maxDisp);

                nn     = cellfun('length', nbrIdx);
                if nnz(nn) > 0
                    rows0 = repelem(find(keepTrack), nn);     % active rows (local)
                    cols0_all = find(keepSpot);
                    cols0 = cols0_all([nbrIdx{:}]');          % spot rows (local)
    
                    % Squared distances (no sqrt)
                    dx  = activeXY(rows0,1) - xyzThis(cols0,1);
                    dy  = activeXY(rows0,2) - xyzThis(cols0,2);
                    d2  = dx.*dx + dy.*dy;

                    % Compact row/col maps
                    [uRows,~,rmap] = unique(rows0);  
                    R = numel(uRows);
                    [uCols,~,cmap] = unique(cols0);  
                    C = numel(uCols);
                    
                    %---- 3) Sparse LAP on rectangular cost with gating ----
                    S = sparse(rmap, cmap, d2, R, C);
                    pairs_sub = matchpairs(S, th2);          % [r c] in compact indices
                    
                    % Back to local active/spot indices, then merge with RNN
                    if ~isempty(pairs_sub)
                        pairs_sub = [uRows(pairs_sub(:,1)), uCols(pairs_sub(:,2))];
                        pairs = [pairs_rnn; pairs_sub];
                    else
                        pairs = pairs_rnn;
                    end
                else
                    pairs = pairs_rnn;
                end
            else
                pairs = pairs_rnn;
            end
        end

        % Validate distance requirement for pairing
        if ~isempty(pairs)
            tr = pairs(:,1);
            sp = pairs(:,2);            

            % Check the distance against the threshold
            idxThis = edges(tIdx):edges(tIdx+1)-1;
            dx = activeXY(tr,1) - x(idxThis(sp));
            dy = activeXY(tr,2) - y(idxThis(sp));
            d2_assigned = dx.^2 + dy.^2;
            
            if any(~all(d2_assigned <= (th2 + 1e-12)))
                ok = d2_assigned <= (th2 + 1e-12); % Links that are less than the distance threshold
                pairs = pairs(ok,:); % drop invalid links
            end
        end



        linkedTracks = false(nActive,1);
        linkedSpots  = false(nNow,1);
    
        % ---- Apply assignments ----
        % After forming final 'pairs' for a frame:
        if ~isempty(pairs)
            tr = pairs(:,1);               % rows in active arrays
            sp = pairs(:,2);               % rows within idxThis
    
            linkedTracks(tr) = true;
            linkedSpots(sp)  = true;
    
            rowIdxVec = idxThis(sp);
            tidVec    = activeID(tr);
            trackID(rowIdxVec) = tidVec;
    
            activeXY(tr,:)  = xyzThis(sp,:);
            activeFrame(tr) = t;
            activeGap(tr)   = 0;
        end
    
        % ---- Start new tracks for unassigned detections ----
        newSpots = find(~linkedSpots);
        if ~isempty(newSpots)
            nNew   = numel(newSpots);
            newIDs = nextID + uint32(0:nNew-1);
            rowIdx = idxThis(newSpots);
    
            trackID(rowIdx) = newIDs;
    
            activeID    = [activeID;   newIDs(:)];
            activeXY    = [activeXY;   xyzThis(newSpots,:)];
            activeFrame = [activeFrame; repmat(t,[nNew 1])];
            activeGap   = [activeGap;  zeros(nNew,1)];
            nextID = nextID + uint32(nNew);
        end
    
        % ---- Increment gaps and terminate old tracks ----
        if nActive>0
            incMask = ~linkedTracks;
            activeGap(incMask) = activeGap(incMask) + 1;
    
            termMask = activeGap > maxGap;
            if any(termMask)
                keep = ~termMask;
                activeID    = activeID(keep);
                activeXY    = activeXY(keep,:);
                activeFrame = activeFrame(keep);
                activeGap   = activeGap(keep);
            end
        end
    end
    
    % -------------------------- Build outputs --------------------------
    % Single, vectorized table write-back
    spotsTbl.TrackID = trackID;
    
    % Group rows by track once
    nTracks = double(max(trackID));
    if nTracks==0
        tracks = struct('id',{},'frames',{},'x',{},'y',{}, ...
                        'quality',{},'Peak',{},'Net_Signal',{}, ...
                        'Integrate_Signal',{},'Integrate_bkgd',{});
        return
    end
    
    have = find(trackID>0);
    rowsByTrack = accumarray(double(trackID(have)), have, [nTracks 1], @(v){v}, {[]});
    
    % Preallocate track struct
    tracks = repmat(struct('id',[], 'frames',[], 'x',[], 'y',[], ...
                           'quality',[], 'Peak',[], 'Net_Signal',[], ...
                           'Integrate_Signal',[], 'Integrate_bkgd',[]), nTracks,1);

    % Fill tracks from arrays (no table access)
    for tid = 1:nTracks
        rows = rowsByTrack{tid};
        if isempty(rows), continue; end
        tracks(tid).id      = tid;
        tracks(tid).frames  = fr(rows).';
        tracks(tid).x       = x(rows).';
        tracks(tid).y       = y(rows).';
        tracks(tid).quality = q(rows).';
        tracks(tid).Peak    = P(rows).';
        tracks(tid).Net_Signal        = NS(rows).';
        tracks(tid).Integrate_Signal  = IS(rows).';
        tracks(tid).Integrate_bkgd    = IB(rows).';
    end
end