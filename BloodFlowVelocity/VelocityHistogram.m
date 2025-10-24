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
    name {mustBeTextScalar, mustBeMember(name, {'artery', 'vein'})}
    n (1, 1) {mustBeInteger, mustBePositive} = 256
end

% Get ToolBox parameters
ToolBox = getGlobalToolBox;
params = ToolBox.getParams;
exportVideos = params.exportVideos;
[~, ~, numFrames] = size(v_video);

% Pre-process data
v_histo = v_video .* mask;
v_min = min(v_histo(mask));
v_max = max(v_histo(mask));

% Set colormap
if strcmp(name, 'artery')
    cmap = ToolBox.Cache.cmapArtery;
else % vein
    cmap = ToolBox.Cache.cmapVein;
end

% Velocity Histogram
T = ToolBox.stride / (1000 * ToolBox.fs); % time between frames
xAx = [0 numFrames * T]; % time axis
yAx = [v_min v_max];

% Initialize histogram matrix
histo = zeros(n, numFrames);
edges = linspace(v_min, v_max, n + 1); % edges for the histcount bins so n+1 numbers

% Initialize figure
fDistrib = figure("Visible", 'on', 'Color', 'w');
fDistrib.Position(3:4) = [600 275];

parfor frameIdx = 1:numFrames
    data = v_histo(:, :, frameIdx);
    histo(:, frameIdx) = histcounts(data(mask), edges); % histcount is faster than histogram or manual for loop counting
end

imagesc(xAx, yAx, histo(:, :));
set(gca, 'PlotBoxAspectRatio', [2.5 1 1])
set(gca, 'YDir', 'normal')
colormap(cmap)
ylabel('Velocity (mm/s)')
xlabel('Time (s)')
histoVideo = frame2im(getframe(fDistrib));
hold on

if exportVideos
    RGB = fill([0 0 xAx(2) xAx(2)], [yAx(1) yAx(2) yAx(2) yAx(1)], [0 0 0], 'EdgeColor', 'none');
    [numX_fig, numY_fig, ~] = size(histoVideo);
    histoVideo = zeros(numX_fig, numY_fig, 3, numFrames);
    gifWriter = GifWriter(sprintf("histogramVelocity%s", name), numFrames);

    for frameIdx = 1:numFrames
        RGB.XData = [frameIdx frameIdx numFrames numFrames] * T;
        f = getframe(fDistrib);
        histoVideo(:, :, :, frameIdx) = frame2im(f);
        gifWriter.write(f, frameIdx);
    end

    gifWriter.generate();
    gifWriter.delete();
end

histoVideo = mat2gray(histoVideo);

exportgraphics(gca, fullfile(ToolBox.path_png, sprintf("%s_histogramVelocity%s.png", ToolBox.folder_name, name)));
exportgraphics(gca, fullfile(ToolBox.path_eps, sprintf("%s_histogramVelocity%s.eps", ToolBox.folder_name, name)));

close(fDistrib)

end
