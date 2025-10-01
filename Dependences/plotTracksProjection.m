function plotTracksProjection(img_stack, tracks, varargin)
%PLOTTRACKSPROJECTION  Show max‑intensity projection with coloured tracks.
%
%   plotTracksProjection(img_stack, tracksStruct)
%
%   INPUTS
%   ------
%   img_stack : H×W×T numeric array (grayscale stack).
%   tracks    : structure array returned by trackSpots() with fields:
%                 id (scalar), x, y, frames   – vectors of equal length.
%
%   OPTIONAL NAME‑VALUE PAIRS
%     'LineWidth' : positive scalar      – line thickness   (default 1.2)
%     'Alpha'     : 0‒1 scalar           – line transparency (default 0.9)
%     'CMap'      : N×3 double colormap  – override colour set (default lines)
%
%   The function displays a figure; it does not return outputs. It assumes
%   the current axes or creates a new one.
% -------------------------------------------------------------------------
% Kevin Scrudders – July 2025
% -------------------------------------------------------------------------

    % Parse inputs
    p = inputParser;
    p.addParameter('LineWidth', 1.2, @(x) isnumeric(x) && isscalar(x) && x>0);
    p.addParameter('Alpha',     0.9, @(x) isnumeric(x) && isscalar(x) && x>=0 && x<=1);
    p.addParameter('CMap',      [],  @(x) isempty(x) || (ismatrix(x) && size(x,2)==3));
    p.parse(varargin{:});
    Lw    = p.Results.LineWidth;
    Alpha = p.Results.Alpha;
    Cmap  = p.Results.CMap;
    edgeColor = [0 0 0];
    
    % Maximum‑intensity projection
    proj = max(img_stack, [], 3);
    proj = squeeze(proj); % ensure 2‑D
    
    % Display projection
    imshow(proj,[min(proj,[],'all') prctile(proj,99.5,'all')]);
    hold on
    
    % Prepare colours
    nTracks = numel(tracks);
    if isempty(Cmap)
        baseCmap = lines(max(nTracks, 7));
    else
        baseCmap = Cmap;
    end
    
    %-------------------- Overlay trajectories + points -------------------
    for k = 1:nTracks
        tk = tracks(k);

        % Skip degenerate or empty tracks
        if ~isfield(tk,'x') || ~isfield(tk,'y') || numel(tk.x) < 1 || numel(tk.x) ~= numel(tk.y)
            continue
        end

        % Color for this track (wrap if fewer rows than tracks)
        c = baseCmap( mod(k-1, size(baseCmap,1)) + 1, : );

        % 1) Scatter every observed (x,y) point for this track
        s = scatter(tk.x, tk.y, 15, c, 'filled');
        % Per-point transparency is well-supported:
        s.MarkerFaceAlpha = Alpha;
        s.MarkerEdgeAlpha = Alpha;

        if ~ischar(edgeColor) && ~isstring(edgeColor)
            s.MarkerEdgeColor = edgeColor;
        else
            s.MarkerEdgeColor = 'none';
        end

        % 2) Connect samples with a line when there are ≥2 points
        if numel(tk.x) >= 2
            % Try RGBA (newer MATLAB), else fall back to opaque RGB
            try
                plot(tk.x, tk.y, '-', 'LineWidth', Lw, 'Color', [c AlphaLine]);
            catch
                plot(tk.x, tk.y, '-', 'LineWidth', Lw, 'Color', c);
            end
        end
    end

    %{
    % Overlay each trajectory
    for k = 1:nTracks
        tk = tracks(k);
        if isempty(tk.x) || numel(tk.x) < 2
            % Skip degenerate tracks
            continue;
        end
        c = baseCmap(mod(k-1, size(baseCmap,1)) + 1, :);
        if Alpha < 1
            c = [c Alpha];          % assume modern MATLAB supports RGBA
        end
        plot(tk.x, tk.y, '-', 'LineWidth', Lw, 'Color', c);
    end
    %}
    
    hold off
    title('Tracks on Maximum‑Intensity Projection');
end
