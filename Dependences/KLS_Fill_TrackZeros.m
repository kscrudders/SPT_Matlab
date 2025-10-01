function out = KLS_Fill_TrackZeros(tracks)
%KLS_Fill_TrackZeros  Linearly interpolate internal gaps in track rows.
%
% INPUT
%   tracks : N x T x C numeric (typically C=2 for x,y). Missing samples may
%            be encoded as NaN, or as (all-zero) vectors across channels.
%
% OUTPUT
%   out : same size as tracks, with INTERNAL gaps (between first and last
%         valid sample per track) linearly interpolated. Leading/trailing
%         gaps are preserved (left as NaN or zeros as in the input).
%
% RULES
%   • A sample at frame t is considered "valid" iff all channels at t are
%     finite and not all zeros simultaneously.
%   • Only frames between the first and last valid samples are candidates
%     for interpolation. Existing valid samples are left unchanged.
%
% NOTES
%   • This is vector-safe (no 1-by-1 indexing pitfalls), robust to multiple
%     disjoint gaps, and to runs of consecutive missing frames.

    out = tracks;                     % preserve dtype/shape
    [N,T,C] = size(tracks);

    % Work track-by-track; per-track work is O(T*C)
    for i = 1:N
        % Xi: T x C view of this track
        Xi = permute(tracks(i,:,:), [2 3 1]);  % (T x C)

        % Missing if any channel is non-finite OR all channels are exactly zero
        nonFinite = any(~isfinite(Xi), 2);     % (T x 1)
        allZero   = all(Xi == 0, 2);           % (T x 1)
        miss      = nonFinite | allZero;       % (T x 1)
        valid     = ~miss;                     % (T x 1)

        vIdx = find(valid);
        if numel(vIdx) < 2
            % Not enough anchors to interpolate — nothing to do
            continue;
        end

        % We only interpolate inside [firstValid, lastValid]
        tFill = vIdx(1):vIdx(end);             % contiguous span
        fillMask = miss(tFill);                % only fill where originally missing

        if ~any(fillMask)
            continue;                          % nothing missing in the interior
        end

        % For each channel, compute linear interpolation over tFill,
        % using only timestamps with valid data as anchors.
        for k = 1:C
            yValid  = Xi(vIdx, k);             % data at valid frames
            % Linear interpolation; no extrapolation needed because tFill
            % lies within [vIdx(1), vIdx(end)].
            yInterp = interp1(vIdx, yValid, tFill, 'linear');

            % Write back only at frames that were missing
            out(i, tFill(fillMask), k) = yInterp(fillMask);
        end
    end
end