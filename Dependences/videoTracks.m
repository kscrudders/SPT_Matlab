function videoTracks(img_stack, tracks, outFile, varargin)
% videoTracks_overlay  Movie with track overlays (format/style like IRCE_gen_Nch_movie).
%
% INPUT
%   img_stack : H x W x T  numeric
%   tracks    : N x T x 2  (N tracks, T frames, 3rd dim: 1=x, 2=y; NaN = missing)
%   outFile   : char/string, e.g. 'tracks.mp4'
%
% NAME–VALUE PAIRS
%   'FPS'        : 20
%   'CLim'       : []                % [min max]; [] = global
%   'Tail'       : Inf                 % # past frames to draw (<= t)
%   'Fade'       : 1                 % fade toggle flag for track tails to decrease in opacity. 1 >= fade, 0 = no fade
%   'MarkerSize' : 6              % scalebar text size
%   'LineWidth'  : 1.5           % track line width
%   'Colors'     : []                % N x 3 RGB; [] = lines(N)
%   'id_specific_Tracks' : []  % List of track numbers that will be colored specifically
%   'Codec'      : 'MPEG-4'     % File format, I suggest staying with mp4
%   'PixelSize'  : 0.157         % µm/px, for optional scalebar
%   'ScaleBar'   : 10             % µm length for scalebar; [] = no bar
%   'Figure'     : []               % reuse axes handle
%   'HideAfterEnd' : 1         % hide track tails after the track ends
%
% KLS style: imshow(...,[0 1]) etc. No external helpers required.
close all

% ---------- options ----------
p = inputParser;
p.addParameter('FPS', 20);
p.addParameter('CLim', []);
p.addParameter('Tail', []);
p.addParameter('Fade', 1);
p.addParameter('MarkerSize', 5);
p.addParameter('LineWidth', 1.5);
p.addParameter('Colors', []);
p.addParameter('id_specific_Tracks',  [], @(x) isnumeric(x) && isvector(x));   % color particular tracks one color
p.addParameter('Codec', 'MPEG-4');
p.addParameter('PixelSize', 0.157);
p.addParameter('ScaleBar', 10);
p.addParameter('Figure', []);
p.addParameter('HideAfterEnd', 1);
p.parse(varargin{:});
o = p.Results;

o.Tail = round(size(img_stack,3)/4);
if o.Fade > 0
    o.Fade = o.Tail;
end

[H,W,T] = size(img_stack);
[N,Tt,C] = size(tracks);
assert(Tt==T && C==2, 'tracks must be N x T x 2');

% colors
if isempty(o.Colors)
    % default: everything blue
    o.Colors = repmat([0 0 1], N, 1);           % pure blue

    % override chosen indices with red
    redIdx = intersect(o.id_specific_Tracks(:), 1:N);    % keep valid indices only
    if ~isempty(redIdx)
        o.Colors(redIdx, :) = repmat([1 0 0], numel(redIdx), 1);
    end
else
    % user‑supplied colour matrix takes precedence over id_specific_Tracks
    assert(size(o.Colors,1) == N && size(o.Colors,2) == 3, ...
        '''Colors'' must be N×3');
end

% intensity limits
if isempty(o.CLim)
    %vmin = double(min(img_stack(:)));
    vmin = double(prctile(img_stack,45,'all')); % works ok with sparse fluorescent data
    vmax = double(prctile(img_stack,99.9,'all')); % works ok with sparse fluorescent data
else
    vmin = o.CLim(1); 
    vmax = o.CLim(2);
end

% figure/axes
if isempty(o.Figure)
    fig = figure('Visible','on','Color','k');
    ax  = axes(fig);
else
    ax  = o.Figure;
    fig = ancestor(ax,'figure');
end

imh = imshow(img_stack(:,:,1), [vmin vmax], 'Parent', ax, ...
             'Border','tight','InitialMagnification',200);
hold(ax,'on'); 
%axis(ax,'image'); 
set(ax,'Visible','off');

% ---------- video writer setup ----------
vidObj = VideoWriter(outFile, o.Codec);
vidObj.Quality = 100;
close(vidObj);

max_t_s = 15;
vidObj.FrameRate = round(min(ceil(T/ max_t_s), 60)); % dynamically set the frame rate to have a video length of 15s maximum
open(vidObj);
drawnow

% --- scalebar (optional) ---
if ~isempty(o.PixelSize) && ~isempty(o.ScaleBar)
    pxLen = round(o.ScaleBar / o.PixelSize);          % pixels
    x0 = W - pxLen - 10;                               % left end
    y0 = H - 10;                                       % line y (near bottom)

    sb = line(ax, [x0 x0+pxLen], [y0 y0], 'Color',[1 1 1], 'LineWidth',3);

    % Text centered horizontally on the bar, placed *below* it
    xMid   = x0 + pxLen/2;
    dyText = 1;                                        % pixels downward (y increases downward)
    text(ax, xMid, y0 + dyText, sprintf('%g µm', o.ScaleBar), ...
        'Color','w','FontSize',10,'FontWeight','bold', ...
        'HorizontalAlignment','center','VerticalAlignment','top');
else
    sb = [];
end


if o.Fade == 0
    % ---------- graphic handles for line plots for tails ----------
    lineH = gobjects(N,1); markH = gobjects(N,1);
    for i = 1:N
        lineH(i) = plot(ax, nan, nan, '-', 'LineWidth', o.LineWidth, ...
                        'Color', o.Colors(i,:));
        markH(i) = plot(ax, nan, nan, 'o', 'MarkerSize', o.MarkerSize, ...
                        'MarkerFaceColor', o.Colors(i,:), ...
                        'MarkerEdgeColor', 'k');
    end
else
    % ---------- graphic handles for patchs instead of line plots for tails ----------
    % this enables opacity
    patchH = gobjects(N,1);   % tails (with alpha gradient)
    markH  = gobjects(N,1);   % current positions
    
    for i = 1:N
        patchH(i) = patch(ax, nan, nan, 'w', ...           % dummy init
            'FaceColor','none', ...
            'EdgeColor',o.Colors(i,:), ...
            'LineWidth',o.LineWidth, ...
            'FaceVertexAlphaData',[], ...
            'AlphaDataMapping','none', ...
            'EdgeAlpha','interp', ...
            'Visible','off');
    
        markH(i) = plot(ax, nan, nan, 'o', ...
            'MarkerSize', o.MarkerSize, ...
            'MarkerFaceColor', o.Colors(i,:), ...
            'MarkerEdgeColor', 'k', ...
            'Visible','off');
    end    
end

% ---------- track end lookup ----------
if o.HideAfterEnd == 1
    % last frame where track i exists (both x & y non-NaN)
    lastValid = zeros(N,1);
    for i = 1:N
        idx = find(~isnan(tracks(i,:,1)) & ~isnan(tracks(i,:,2)), 1, 'last');
        if isempty(idx), idx = 0; end
        lastValid(i) = idx;
    end
else
    lastValid = ones(N,1);
    lastValid =  lastValid .* T;
end

% ---------- track tail fade lookup ----------
Ltail = o.Tail;                          % number of frames we keep
if o.Fade <= 0
    fadeLUT = ones(Ltail,1);
else
    % linear example; replace with exp(-(0:Ltail-1)'/tau) for exponential
    fadeLUT = max(0, 1 - (0:Ltail-1)'/o.Fade);
end
fadeLUT = flip(fadeLUT);

lastImg = [];  % initialize prior image slice (outside the loop)

% ---------- main loop ----------
for t = 1:T
    currentImg = img_stack(:,:,t);

    % Skip frame if identical to previous
    if ~isempty(lastImg) && isequal(currentImg, lastImg)
        drawnow limitrate;
        writeVideo(vidObj, getframe(fig));
        continue
    end

    lastImg = currentImg;

    % ---------- update image ----------
    imh.CData = currentImg;

    % ---------- update track tail range ----------
    % tail frame range
    if isinf(o.Tail), fStart = 1; else, fStart = max(1, t - o.Tail); end
    fRange = fStart:t;

    Xt = tracks(:,fRange,1);  % N x L
    Yt = tracks(:,fRange,2);
    Xc = tracks(:,t,1);
    Yc = tracks(:,t,2);

    % ---------- per-track visualization block ----------
    if o.Fade == 0
        % ---------- line plots for tails (no opacity) ----------
        for i = 1:N
            if t > lastValid(i) || lastValid(i) == 0
                % track is over (or never existed) -> hide
                set(lineH(i),'Visible','off');
                set(markH(i),'Visible','off');
                continue
            end
        
            % Tail range limited by lastValid(i)
            fStart = max(1, t - o.Tail);
            fRange_i = fStart : t;
            fRange_i(fRange_i > lastValid(i)) = [];
        
            xi = tracks(i,fRange_i,1);
            yi = tracks(i,fRange_i,2);
            good = ~isnan(xi) & ~isnan(yi);
        
            set(lineH(i),'XData',xi(good),'YData',yi(good), ...
                          'Color',o.Colors(i,:),'Visible','on');
        
            % current marker
            Xc = tracks(i,t,1);  Yc = tracks(i,t,2);
            if ~isnan(Xc) && ~isnan(Yc)
                set(markH(i),'XData',Xc,'YData',Yc,'Visible','on');
            else
                set(markH(i),'Visible','off');
            end
        end
    else
        % ---------- patch lines for tails (with opacity) ----------
        for i = 1:N
            if t > lastValid(i) || lastValid(i) == 0
                set(patchH(i),'Visible','off');
                set(markH(i),'Visible','off');
                continue
            end
        
            % Tail indices limited by last valid frame
            fStart   = max(1, t - Ltail + 1);
            fRange_i = fStart : t;
            fRange_i(fRange_i > lastValid(i)) = [];
        
            xi = tracks(i,fRange_i,1);
            yi = tracks(i,fRange_i,2);
            good = ~isnan(xi) & ~isnan(yi);
            xi = xi(good);  yi = yi(good);
            if isempty(xi)
                set(patchH(i),'Visible','off');
                set(markH(i),'Visible','off');
                continue
            end
        
            % Alpha vector (same length as vertices)
            nSeg = numel(xi);                     % vertices count
            alphaVec = fadeLUT(end-nSeg+1:end);   % youngest last element = 1
            % Ensure column vector
            alphaVec = alphaVec(:);
        
            % Break the line so patch doesn't wrap
            Xv = [xi(:); NaN];
            Yv = [yi(:); NaN];
            Av = [alphaVec(:); NaN];   % NaN (or 0) is fine as a terminator

            % Update patch
            set(patchH(i), 'XData', Xv, 'YData', Yv, ...
                           'FaceVertexAlphaData', Av, ...
                           'Visible','on');
            %{
            set(patchH(i), 'XData', xi, 'YData', yi, ...
                           'FaceVertexAlphaData', alphaVec, ...
                           'Visible','on');
            %}
            % Marker at current frame
            Xc = tracks(i,t,1);  Yc = tracks(i,t,2);
            if ~isnan(Xc) && ~isnan(Yc)
                set(markH(i),'XData',Xc,'YData',Yc,'Visible','on');
            else
                set(markH(i),'Visible','off');
            end
        end
    end


    % ---------- update video writer ----------
    %drawnow limitrate nocallbacks;
    drawnow limitrate;
    writeVideo(vidObj, getframe(fig));
end

close(vidObj);
if isempty(o.Figure), close(fig); end
end
