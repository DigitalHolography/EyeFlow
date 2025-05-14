function histoVideo = VelocityHistogram(v_video, mask, name, n)
% VelocityHistogram.m: Creates an Histogram of the velocities from a video
% and passes the figure as a mat
%
% Inputs:
%   v_video     : 3D array of velocities (x,y,frame)
%   mask        : Logical mask to select velocities
%   name        : 'Arteries' or 'Veins' for naming/coloring
%   n           : Resolution of the histogram (default: 256)
%
% Outputs:
%   histoVideo  : Matrix of the histogram figure

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

% Pre-process data
v_histo = v_video .* mask;
validPixels = mask ~= 0;
v_min = min(v_histo(validPixels));
v_max = max(v_histo(validPixels));

% Set colormap
if strcmp(name, 'Artery')
    cmap = ToolBox.cmapArtery;
else % vein
    cmap = ToolBox.cmapVein;
end

% Initialize figure
fig = figure("Visible", 'off', 'Color', 'w');
fig.Position(3:4) = [600 275];
ax = gca;

% Set up axes and appearance
xAx = [0, size(v_video, 3) * ToolBox.stride / (1000 * ToolBox.fs)];
yAx = [v_min, v_max];
set(ax, 'YDir', 'normal', 'PlotBoxAspectRatio', [2.5 1 1]);
fontsize(ax, 14, "points");
colormap(ax, cmap);
ylabel(ax, 'Velocity (mm.s^{-1})');
xlabel(ax, 'Time (s)');
title(ax, sprintf("Velocity distribution in %s", name));

% Initialize histogram matrix
histo = zeros(n, size(v_video, 3));
edges = linspace(v_min, v_max, n + 1);

% Pre-compute bin indices for all frames
if exportVideos
    % For video export: process frame by frame
    histoVideo = zeros(fig.Position(4), fig.Position(3), 3, size(v_video, 3));

    for frameIdx = 1:size(v_video, 3)
        frameData = v_histo(:, :, frameIdx);
        validData = frameData(validPixels(:, :, min(frameIdx, size(validPixels, 3))));

        % Vectorized histogram calculation
        [counts, ~] = histcounts(validData, edges);
        histo(:, frameIdx) = counts';

        % Update display
        imagesc(ax, xAx, yAx, histo);
        drawnow;

        set(ax, 'YDir', 'normal', 'PlotBoxAspectRatio', [2.5 1 1]);
        fontsize(ax, 14, "points");
        colormap(ax, cmap);
        ylabel(ax, 'Velocity (mm.s^{-1})');
        xlabel(ax, 'Time (s)');
        title(ax, sprintf("Velocity distribution in %s", name));

        % Capture frame
        histoVideo(:, :, :, frameIdx) = frame2im(getframe(fig));
    end

    writeGifOnDisc(mat2gray(histoVideo), sprintf("histogramVelocity%s", name), "ToolBox", ToolBox);
else
    % For single image: process all frames at once
    for frameIdx = 1:size(v_video, 3)
        frameData = v_histo(:, :, frameIdx);
        validData = frameData(validPixels(:, :, min(frameIdx, size(validPixels, 3))));

        % Vectorized histogram calculation
        [counts, ~] = histcounts(validData, edges);
        histo(:, frameIdx) = counts';
    end

    % Display final histogram
    imagesc(ax, xAx, yAx, histo);
    drawnow;

    set(ax, 'YDir', 'normal', 'PlotBoxAspectRatio', [2.5 1 1]);
    fontsize(ax, 14, "points");
    colormap(ax, cmap);
    ylabel(ax, 'Velocity (mm.s^{-1})');
    xlabel(ax, 'Time (s)');
    title(ax, sprintf("Velocity distribution in %s", name));

    % Capture final image
    histoVideo = frame2im(getframe(fig));
end

histoVideo = mat2gray(histoVideo);

% Export graphics
outputDir = fullfile(ToolBox.path_png, 'bloodFlowVelocity');

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

exportgraphics(ax, fullfile(outputDir, sprintf("%s_histogramVelocity%s.png", ToolBox.main_foldername, name)));

outputDir = fullfile(ToolBox.path_eps, 'bloodFlowVelocity');

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

exportgraphics(ax, fullfile(outputDir, sprintf("%s_histogramVelocity%s.eps", ToolBox.main_foldername, name)));

fprintf("- Velocity Histogram %s Timing: %.2fs\n", name, toc(tVelocityVideo));
end
