function [videoPlotMean, videoPlotFrames] = videoPlot(video, mask, xy_barycenter, v_mean_RGB, v_video_RGB, opt)

arguments
    video
    mask
    xy_barycenter
    v_mean_RGB
    v_video_RGB
    opt.etiquettes_locs = []
    opt.etiquettes_values = []
    opt.Visible = "off"
end

% Get global toolbox settings
ToolBox = getGlobalToolBox;
params = ToolBox.getParams;
exportVideos = params.exportVideos;

% Rescale video and get dimensions
[numX, numY, numFrames] = size(video);
fontsize = round(numX / 75);
video = rescale(video);

% Precompute constants
x_c = xy_barycenter(1);
y_c = xy_barycenter(2);
r2 = params.json.SizeOfField.BigRadiusRatio;
r1 = params.json.SizeOfField.SmallRadiusRatio;
w = 0.005;

% Precompute circles
circle1 = diskMask(numX, numY, r1, r1 + w, center = [x_c / numX, y_c / numY]);
circle2 = diskMask(numX, numY, r2 - w, r2, center = [x_c / numX, y_c / numY]);
maskCircles = circle1 | circle2;

% Initialize video plot
locs = opt.etiquettes_locs;
values = opt.etiquettes_values;

% Initialize signal plot with optimized properties
figure('Name', 'Signal Plot', 'Color', 'w', ...
    'Visible', opt.Visible, ...
    'Position', [200 200 600 600]);

image = rescale(mean(video, 3));
image_RGB = v_mean_RGB .* mask + image .* ~mask;
image_RGB = image_RGB .* ~(maskCircles & ~mask) + maskCircles .* ~mask;
imshow(image_RGB);

axis image
axis off

t = cell(1, size(locs, 1));
new_x = cell(1, size(locs, 1));
new_y = cell(1, size(locs, 1));

meanQ = mean(values, 2);

if ~isempty(locs)

    for etIdx = 1:size(locs, 1)
        new_x{etIdx} = locs(etIdx, 1);
        new_y{etIdx} = locs(etIdx, 2);

        % Add the text
        t{etIdx} = text(new_x{etIdx}, new_y{etIdx}, sprintf('%0.1f', meanQ(etIdx)), ...
            "FontWeight", "bold", ...
            "FontSize", fontsize, ...
            "Color", "white", ...
            "BackgroundColor", "black");

    end

end

videoPlotMean = frame2im(getframe(gca));

% Preallocate video data array
videoPlotFrames = zeros([numX, numY, 3, numFrames], 'single');

if exportVideos

    % Generate video frames
    parfor frameIdx = 1:numFrames
        figure('Name', 'Signal Plot', 'Color', 'w', ...
            'Visible', opt.Visible, ...
            'Position', [200 200 600 600]);

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

        frame = frame2im(getframe(gca));

        if size(videoPlotFrames(:, :, :, frameIdx), 1) == size(frame, 1)
            videoPlotFrames(:, :, :, frameIdx) = frame;
        else
            videoPlotFrames(:, :, :, frameIdx) = imresize(frame, [numX, numY]);
        end

    end

end

end
