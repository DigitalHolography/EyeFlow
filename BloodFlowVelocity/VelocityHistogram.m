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

if strcmp(name, 'Artery')
    cmap = ToolBox.cmapArtery;
elseif strcmp(name, 'Vein')
    cmap = ToolBox.cmapVein;
else
    cmap = cmapLAB(256, [0 0 0], 0, [0 1 0], 1/3, [1/2 1 1/2], 2/3, [1 1 1], 1);
end

yAx = [v_min v_max];

%% Velocity Histogram

X = linspace(v_min, v_max, n);
xAx = [0 numFrames * ToolBox.stride / (1000 * ToolBox.fs)];
histo = zeros(n, numFrames);
D = (v_max - v_min) / (n);
edges = linspace(v_min-D/2,v_max+D/2,n+1); % edges for the histcount bins so n+1 numbers


fDistrib = figure("Visible", 'on', 'Color', 'w');
fDistrib.Position(3:4) = [600 275];

h_imagesc = imagesc(xAx, yAx, histo(:, :));
set(gca, 'YDir', 'normal')
ylabel('Velocity (mm.s^{-1})')
xlabel('Time (s)')
title(sprintf("velocity distribution in %s", name))
set(gca, 'PlotBoxAspectRatio', [2.5 1 1])
colormap(cmap)
[numX_fig, numY_fig] = deal(fDistrib.Position(4),fDistrib.Position(3));


if exportVideos
    histoVideo = zeros(numX_fig, numY_fig, 3, numFrames);
    gifWriter = GifWriter(sprintf("histogramVelocity%s", name), numFrames);

    for frameIdx = 1:numFrames
        data = v_histo(:,:,frameIdx);
        histo(:, frameIdx) = histcounts(data(mask),edges); % histcount is faster than histogram or manual for loop counting
        
        h_imagesc.CData(:, frameIdx) = histo(:, frameIdx);

        f = getframe(fDistrib);
        histoVideo(:, :, :, frameIdx) = frame2im(f);
        gifWriter.write(f, frameIdx);
    end

    gifWriter.generate();
    gifWriter.delete();

    exportgraphics(gca, fullfile(ToolBox.path_png, 'bloodFlowVelocity', sprintf("%s_histogramVelocity%s.png", ToolBox.main_foldername, name)))
    exportgraphics(gca, fullfile(ToolBox.path_eps, 'bloodFlowVelocity', sprintf("%s_histogramVelocity%s.eps", ToolBox.main_foldername, name)))
else

    for frameIdx = 1:numFrames
        data = v_histo(:,:,frameIdx);
        histo(:, frameIdx) = histcounts(data(mask),edges); % histcount is faster than histogram or manual for loop counting
    end

    h_imagesc.CData(:, :) = histo(:, :);

    exportgraphics(gca, fullfile(ToolBox.path_png, 'bloodFlowVelocity', sprintf("%s_histogramVelocity%s.png", ToolBox.main_foldername, name)))
    exportgraphics(gca, fullfile(ToolBox.path_eps, 'bloodFlowVelocity', sprintf("%s_histogramVelocity%s.eps", ToolBox.main_foldername, name)))
end

fprintf("- Velocity Histogram %s Timing : %ds\n", name, round(toc(tVelocityVideo)))

end
