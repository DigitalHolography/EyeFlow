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
    v_video {mustBeNumeric}
    mask {mustBeNumericOrLogical}
    name {mustBeTextScalar, mustBeMember(name, {'Artery', 'Vein'})}
    n (1, 1) {mustBeInteger, mustBePositive} = 256
end

tVelocityVideo = tic;

% Get toolbox parameters
ToolBox = getGlobalToolBox;
params = ToolBox.getParams;
exportVideos = params.exportVideos;
[~, ~, numFrames] = size(v_video);

% Pre-process data
v_histo = v_video .* mask;
v_min = min(v_histo(mask));
v_max = max(v_histo(mask));

% Set colormap
if strcmp(name, 'Artery')
    cmap = ToolBox.cmapArtery;
else % vein
    cmap = ToolBox.cmapVein;
end

% Velocity Histogram
xAx = [0 numFrames * ToolBox.stride / (1000 * ToolBox.fs)];
yAx = [v_min v_max];

% Initialize histogram matrix
histo = zeros(n, numFrames);
edges = linspace(v_min, v_max, n + 1); % edges for the histcount bins so n+1 numbers

% Initialize figure
fDistrib = figure("Visible", 'on', 'Color', 'w');
fDistrib.Position(3:4) = [600 275];
h_imagesc = imagesc(xAx, yAx, histo(:, :));
set(gca, 'PlotBoxAspectRatio', [2.5 1 1])
set(gca, 'YDir', 'normal')
colormap(cmap)
ylabel('Velocity (mm/s)')
xlabel('Time (s)')
title(sprintf("velocity distribution in %s", name))

% Pre-compute bin indices for all frames
[numX_fig, numY_fig] = deal(fDistrib.Position(4), fDistrib.Position(3));

if exportVideos
    histoVideo = zeros(numX_fig, numY_fig, 3, numFrames);
    gifWriter = GifWriter(sprintf("histogramVelocity%s", name), numFrames);

    for frameIdx = 1:numFrames
        data = v_histo(:, :, frameIdx);
        histo(:, frameIdx) = histcounts(data(mask), edges); % histcount is faster than histogram or manual for loop counting

        h_imagesc.CData(:, frameIdx) = histo(:, frameIdx);

        f = getframe(fDistrib);
        histoVideo(:, :, :, frameIdx) = frame2im(f);
        gifWriter.write(f, frameIdx);
    end

    gifWriter.generate();
    gifWriter.delete();

else

    for frameIdx = 1:numFrames
        data = v_histo(:, :, frameIdx);
        histo(:, frameIdx) = histcounts(data(mask), edges); % histcount is faster than histogram or manual for loop counting
    end

    h_imagesc.CData(:, :) = histo(:, :);
    f = getframe(fDistrib);
    histoVideo = frame2im(f);

end

histoVideo = mat2gray(histoVideo);

% Export graphics
outputDir = fullfile(ToolBox.path_png, 'bloodFlowVelocity');

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

exportgraphics(gca, fullfile(outputDir, sprintf("%s_histogramVelocity%s.png", ToolBox.main_foldername, name)));

outputDir = fullfile(ToolBox.path_eps, 'bloodFlowVelocity');

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

exportgraphics(gca, fullfile(outputDir, sprintf("%s_histogramVelocity%s.eps", ToolBox.main_foldername, name)));

fprintf("- Velocity Histogram %s Timing : %ds\n", name, round(toc(tVelocityVideo)))

end
