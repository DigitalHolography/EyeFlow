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
end

% Get global toolbox settings
ToolBox = getGlobalToolBox;
params = ToolBox.getParams;
exportVideos = params.exportVideos;

% Plot the video
[videoPlotMean, videoPlotFrames] = videoPlot(video, mask, xy_barycenter, v_mean_RGB, v_video_RGB, ...
    'etiquettes_locs', opt.etiquettes_locs, 'etiquettes_values', opt.etiquettes_values);

% Plot the signal graph
[signalPlotMean, signalPlotFrames] = signalPlot(signal, stdsignal);

% Combines the signal graph with the image video
[~, numY_fig, ~, ~] = size(signalPlotFrames);
combinedLast = cat(1, mat2gray(imresize(videoPlotMean, [numY_fig numY_fig])), mat2gray(signalPlotMean));

% Save final frames as PNGs
imwrite(mat2gray(signalPlotMean), fullfile(ToolBox.path_png, 'local', sprintf("%s_%s_plot.png", ToolBox.folder_name, dirname)));
imwrite(mat2gray(videoPlotMean), fullfile(ToolBox.path_png, 'local', sprintf("%s_%s.png", ToolBox.folder_name, dirname)));
imwrite(combinedLast, fullfile(ToolBox.path_png, 'local', sprintf("%s_%s_combined.png", ToolBox.folder_name, dirname)));

% Save as GIF if not skipping frames
if exportVideos
    % Combine video and signal frames
    videoInterp = imresize(mat2gray(videoPlotFrames), [size(signalPlotFrames, 2), size(signalPlotFrames, 2)]);
    videoInterp = max(0, min(videoInterp, 1)); % Ensure values are within [0, 1]
    combinedFrames = cat(1, videoInterp, mat2gray(signalPlotFrames));
    writeGifOnDisc(mat2gray(videoPlotFrames), sprintf("%s", dirname), 0.04);
    writeGifOnDisc(mat2gray(signalPlotFrames), sprintf("%s_plot", dirname), 0.04);
    writeGifOnDisc(combinedFrames, sprintf("%s_combined", dirname), 0.04);
end

end
