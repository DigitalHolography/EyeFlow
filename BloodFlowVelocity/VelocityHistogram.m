function histoVideo = VelocityHistogram(v_video, mask, name, n)
% VelocityHistogram.m: Creates an Histogram of the velocities from a video
% and passes the figure as a mat
%
% Inputs:
%   v_video     :   video of the velocities
%   mask        :   a logical mask to locate the velocities
%   name        :   'Arteries' or 'Veins' to change the name of the figs
%               and the color maps
%   n           :   Resolution of the Histrogram
%
% Outputs:
%   histoVideo  :   Mat of the histogram figure

arguments
    v_video
    mask
    name
    n = 256
end

tVelocityVideo = tic;

ToolBox = getGlobalToolBox;
params = ToolBox.getParams;
exportVideos = params.exportVideos;

[numX, numY, numFrames] = size(v_video);

v_histo = v_video .* mask;
v_min = min(v_histo, [], 'all');
v_max = max(v_histo, [], 'all');

if strcmp(name, 'artery')
    cmap = ToolBox.cmapArtery;
elseif strcmp(name, 'vein')
    cmap = ToolBox.cmapVein;
else
    cmap = cmapLAB(256, [0 0 0], 0, [0 1 0], 1/3, [1/2 1 1/2], 2/3, [1 1 1], 1);
end

X = linspace(v_min, v_max, n);
xAx = [0 numFrames * ToolBox.stride / (1000 * ToolBox.fs)];
yAx = [v_min v_max];
histo = zeros(n, numFrames);
D = (v_max - v_min) / (n - 1);

fig = figure("Visible", 'off', 'Color', 'w');
fig.Position(3:4) = [600 275];

indexMin = find(X == v_min);
indexMax = find(X == v_max);
imagesc(xAx, yAx, histo(indexMin:indexMax, :))
set(gca, 'YDir', 'normal')
set(gca, 'PlotBoxAspectRatio', [2.5 1 1])
fontsize(gca, 14, "points");
colormap(cmap)
f = getframe(gcf);
[numX_fig, numY_fig, ~] = size(f.cdata);

if exportVideos
    histoVideo = zeros(numX_fig, numY_fig, 3, numFrames);

    for frameIdx = 1:numFrames

        for xx = 1:numX

            for yy = 1:numY

                if mask(xx, yy) ~= 0
                    i = find(and(X >= v_histo(xx, yy, frameIdx), X < v_histo(xx, yy, frameIdx) + D));
                    histo(i, frameIdx) = histo(i, frameIdx) + 1;
                end

            end

        end

        imagesc(xAx, yAx, histo(indexMin:indexMax, :))
        set(gca, 'YDir', 'normal')
        ylabel('Velocity (mm.s^{-1})')
        xlabel('Time (s)')
        title(sprintf("velocity distribution in %s", name))
        set(gca, 'PlotBoxAspectRatio', [2.5 1 1])
        histoVideo(:, :, :, frameIdx) = rescale(frame2im(getframe(gcf)));
    end

    writeGifOnDisc(histoVideo, sprintf("histogramVelocity%s", name), "ToolBox", ToolBox);

    exportgraphics(gca, fullfile(ToolBox.path_png, 'bloodFlowVelocity', sprintf("%s_histogramVelocity%s.png", ToolBox.main_foldername, name)))
    exportgraphics(gca, fullfile(ToolBox.path_eps, 'bloodFlowVelocity', sprintf("%s_histogramVelocity%s.eps", ToolBox.main_foldername, name)))
else

    for frameIdx = 1:numFrames

        for xx = 1:numX

            for yy = 1:numY

                if mask(xx, yy) ~= 0
                    i = find(and(X >= v_histo(xx, yy, frameIdx), X < v_histo(xx, yy, frameIdx) + D));
                    histo(i, frameIdx) = histo(i, frameIdx) + 1;
                end

            end

        end

    end

    figure(fig, 'Visible', 'off')
    imagesc(xAx, yAx, histo(indexMin:indexMax, :))
    set(gca, 'YDir', 'normal')
    ylabel('Velocity (mm.s^{-1})')
    xlabel('Time (s)')
    title(sprintf("velocity distribution in %s", name))
    set(gca, 'PlotBoxAspectRatio', [2.5 1 1])
    fontsize(gca, 14, "points");
    f = getframe(gcf);
    histoVideo = f.cdata;
    exportgraphics(gca, fullfile(ToolBox.path_png, 'bloodFlowVelocity', sprintf("%s_histogramVelocity%s.png", ToolBox.main_foldername, name)))
    exportgraphics(gca, fullfile(ToolBox.path_eps, 'bloodFlowVelocity', sprintf("%s_histogramVelocity%s.eps", ToolBox.main_foldername, name)))

end

fprintf("- Velocity Histogram %s Timing : %ds\n", name, round(toc(tVelocityVideo)))

end
