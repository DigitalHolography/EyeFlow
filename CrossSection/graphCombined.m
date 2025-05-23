function graphCombined(video, v_video_RGB, v_mean_RGB, mask, signal, stdsignal, xy_barycenter, dirname, opt)
% Combines a video and a signal into an animated GIF, overlaying them frame by frame.
%
% Inputs:
%   Videofield: Video to be displayed (grayscale).
%   mask: Mask to overlay in red.
%   signal: Signal data to plot.
%   stdsignal: Standard deviation of the signal.
%   xy_barycenter: Barycenter coordinates [x, y].
%   dirname: Directory name for saving outputs.
%   opt: Optional parameters (see below).
%
% Optional Parameters (opt):
%   etiquettes_locs: Locations of tags (empty if not needed).
%   etiquettes_values: Values for tags (empty if not needed).
%   skip: Skip frames (logical, default false).
%   Color: Color for mask overlay.
%   Visible: Visibility of figures (logical, default false).

arguments
    video
    v_video_RGB
    v_mean_RGB
    mask
    signal
    stdsignal
    xy_barycenter
    dirname
    opt.etiquettes_locs = []
    opt.etiquettes_values = []
    opt.skip logical = false
    opt.Visible logical = false
end

% Get global toolbox settings
ToolBox = getGlobalToolBox;
params = ToolBox.getParams;
exportVideos = params.exportVideos;

% Rescale video and get dimensions
[numX, numY, numFrames] = size(video);
T = ToolBox.stride / ToolBox.fs / 1000;
fullTime = linspace(0, numFrames * T, numFrames);
fontsize = round(numX / 75);
video = rescale(video);

% Precompute constants
x_c = xy_barycenter(1);
y_c = xy_barycenter(2);
r2 = params.json.SizeOfField.BigRadiusRatio;
r1 = params.json.SizeOfField.SmallRadiusRatio;
w = 0.005;
Color_std = [0.7, 0.7, 0.7];

% Precompute circles
circle1 = diskMask(numX, numY, r1, r1 + w, center = [x_c / numX, y_c / numY]);
circle2 = diskMask(numX, numY, r2 - w, r2, center = [x_c / numX, y_c / numY]);
maskCircles = circle1 | circle2;

% Set y-axis limits for the signal plot
mean_signal = mean(signal);

% Determine starting frame based on skip option
if opt.skip
    startingFrame = numFrames;
else
    startingFrame = 1;
end

% Initialize video plot
locs = opt.etiquettes_locs;
values = opt.etiquettes_values;

videoPlot = figure(410);
videoPlot.Position = [200 200 600 600];
videoPlot.Visible = opt.Visible;

image = rescale(mean(video, 3));
image_RGB = v_mean_RGB .* mask + image .* ~mask;
image_RGB = image_RGB .* ~(maskCircles & ~mask) + maskCircles .* ~mask;
imshow(image_RGB);

axis image
axis off

t = cell(1, size(locs, 1));
new_x = cell(1, size(locs, 1));
new_y = cell(1, size(locs, 1));

if ~isempty(locs)

    for etIdx = 1:size(locs, 1)
        new_x{etIdx} = locs(etIdx, 1);
        new_y{etIdx} = locs(etIdx, 2);

        % Add the text
        t{etIdx} = text(new_x{etIdx}, new_y{etIdx}, sprintf('%0.1f', mean(values(etIdx), 2)), ...
            "FontWeight", "bold", ...
            "FontSize", fontsize, ...
            "Color", "white", ...
            "BackgroundColor", "black");

    end

end

frame = frame2im(getframe(gca));

if exportVideos
    % Preallocate video data array
    videoPlotFrames = zeros([numX, numY, 3, numFrames], 'single');

    % Generate video frames
    parfor frameIdx = startingFrame:numFrames
        image_RGB = v_video_RGB(:, :, :, frameIdx) .* mask + video(:, :, frameIdx) .* ~mask;
        image_RGB = image_RGB .* ~(maskCircles & ~mask) + maskCircles .* ~mask;
        imshow(image_RGB);

        if ~isempty(locs)

            for etIdx = 1:size(locs, 1)

                % Add the text
                text(new_x{etIdx}, new_y{etIdx}, sprintf('%0.1f', values(etIdx, frameIdx)), ...
                    "FontWeight", "bold", ...
                    "FontSize", fontsize, ...
                    "Color", "white", ...
                    "BackgroundColor", "black");

            end

        end

        videoPlotFrames(:, :, :, frameIdx) = frame2im(getframe(gca));
    end

end

% Initialize signal plot
signalPlot = figure('Name', 'Signal Plot', "Color", 'w');
signalPlot.Visible = opt.Visible;
signalPlot.Position = [200 200 600 300];

% Generate signal plot frames
curve1 = signal + stdsignal;
curve2 = signal - stdsignal;
tmp_fullTime = [fullTime, fliplr(fullTime)];
inBetween = [curve1, fliplr(curve2)];
fill(tmp_fullTime, inBetween, Color_std);

ylabel('Volume Rate (µL/min)')
xlabel('Time (s)')
title(sprintf("Average Blood Volume Rate : %.0f %s", round(mean_signal), 'µL/min'))
hold on;
plot(fullTime, curve1, "Color", Color_std, 'LineWidth', 2);
plot(fullTime, curve2, "Color", Color_std, 'LineWidth', 2);
plot(fullTime, signal, '-k', 'LineWidth', 2);
yline(mean_signal, '--k', 'LineWidth', 2)

yAx = [min(-1, min(signal)), max(signal) * 1.07];
axss = [0, numFrames * T, min(-1, min(signal)), max(signal) * 1.07];
axis(axss);
box on
set(gca, 'Linewidth', 2)
pbaspect([2.5, 1, 1])

signalPlotFrames = frame2im(getframe(signalPlot));
cache = fill([0 0 numFrames numFrames] * T, [yAx(1), yAx(2), yAx(2), yAx(1)], 'w', "EdgeColor", 'none');
hold off;

if exportVideos
    % Preallocate signal plot data array
    signalPlotFrames = zeros([size(signalPlotFrames), numFrames], 'single');

    for frameIdx = startingFrame:numFrames
        % Update signal plot for each frame
        cache.XData = [frameIdx * T, frameIdx * T, numFrames * T, numFrames * T];
        axis(axss);

        set(gca, 'Linewidth', 2)
        pbaspect([2.5, 1, 1])

        signalPlotFrames(:, :, :, frameIdx) = frame2im(getframe(signalPlot));
    end

end

hold off;

[~, numY_fig, ~, ~] = size(signalPlotFrames);
combinedLast = cat(1, mat2gray(imresize(frame, [numY_fig numY_fig])), mat2gray(signalPlotFrames(:, :, :, end)));

% Save final frames as PNGs
imwrite(mat2gray(signalPlotFrames(:, :, :, end)), fullfile(ToolBox.path_png, 'crossSectionsAnalysis', sprintf("%s_%s_plot.png", ToolBox.folder_name, dirname)));
imwrite(combinedLast, fullfile(ToolBox.path_png, 'crossSectionsAnalysis', sprintf("%s_%s_combined.png", ToolBox.folder_name, dirname)));

% Save as GIF if not skipping frames
if exportVideos
    % Combine video and signal frames
    videoInterp = imresize(mat2gray(videoPlotFrames), [size(signalPlotFrames, 2), size(signalPlotFrames, 2)]);
    videoInterp = max(0, min(videoInterp, 1)); % Ensure values are within [0, 1]
    combinedFrames = cat(1, videoInterp, mat2gray(signalPlotFrames));
    writeGifOnDisc(mat2gray(signalPlotFrames), sprintf("%s_plot", dirname), 0.04);
    writeGifOnDisc(combinedFrames, sprintf("%s_combined", dirname), 0.04);
end

% Close figures
close(videoPlot, signalPlot);
end
