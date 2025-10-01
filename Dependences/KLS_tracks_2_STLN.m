function [STLN_tracks, STLN_intensity] = KLS_tracks_2_STLN(img_stack, tracks)
    N = length(tracks);
    t = size(img_stack,3);

    STLN_tracks = NaN([N t 2]); % empty Structure
    STLN_intensity = NaN([N t 2]); % empty Structure
    
    for i = 1:length(tracks)
        STLN_tracks(tracks(i).id,tracks(i).frames,1) = tracks(i).x; % x
        STLN_tracks(tracks(i).id,tracks(i).frames,2) = tracks(i).y; % y
    
        STLN_intensity(tracks(i).id,tracks(i).frames,1) = tracks(i).Peak; % peak of guassian fit
    end
end