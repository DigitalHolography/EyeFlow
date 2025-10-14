function [v_video_RGB, v_mean_RGB] = flowMap(v_video, maskSection, maskArtery, maskVein, M0_ff, xy_barycenter, ToolBox)

params = ToolBox.getParams;
veinsAnalysis = params.veins_analysis;
exportVideos = params.exportVideos;

% Rescale once
M0_ff = rescale(M0_ff);
M0_ff_image = rescale(mean(M0_ff, 3));
[numX, numY, numFrames] = size(v_video);

% Precompute constants
x_c = xy_barycenter(1) / numX;
y_c = xy_barycenter(2) / numY;
r1 = params.json.SizeOfField.BigRadiusRatio;
r2 = params.json.SizeOfField.SmallRadiusRatio;
w = 0.01;

% Precompute circles
circle1 = diskMask(numX, numY, r1, r1 + w, center = [x_c, y_c]);
circle2 = diskMask(numX, numY, r2 - w, r2, center = [x_c, y_c]);
circles = circle1 | circle2;

% Precompute colormaps
cmapArtery = ToolBox.Cache.cmapArtery;

if veinsAnalysis
    cmapVein = ToolBox.Cache.cmapVein;
    cmapAV = ToolBox.Cache.cmapAV;
end

% Precompute v_rescaled
v_max = max(v_video(maskSection));
v_min = 0;
v_rescaled = abs(v_video);
v_rescaled = (v_rescaled - v_min) / v_max;

% Precompute v_mean and v_mean_rescaled
v_mean = squeeze(mean(v_video, 3));
v_mean_rescaled = squeeze(mean(v_rescaled, 3));
v_video_RGB = zeros(numX, numY, 3, numFrames, 'single');

if veinsAnalysis
    % Velocity Images
    velocityIm(v_mean, maskArtery & maskSection, cmapArtery, 'v_Artery', colorbarOn = true);
    velocityIm(v_mean, maskVein & maskSection, cmapVein, 'v_Vein', colorbarOn = true);
    velocityIm(v_mean, (maskArtery | maskVein) & maskSection, turbo, 'v_Vessel', colorbarOn = true);

    % Precompute masks
    maskAV = maskArtery & maskVein;
    maskArtery = maskArtery & ~maskAV;
    maskVein = maskVein & ~maskAV;
else

    % Velocity Images
    velocityIm(v_mean, maskArtery, cmapArtery, 'v_Artery', colorbarOn = true);

end

if veinsAnalysis
    % velocity Colorbars of the gif
    velocityColorbar(cmapArtery, v_min, v_max, 'v_Artery');
    velocityColorbar(cmapVein, v_min, v_max, 'v_Vein');

    % Average Image
    v_mean_RGB = flowMapImg(v_mean_rescaled, {maskArtery, maskVein, maskAV}, {cmapArtery, cmapVein, cmapAV}, background = M0_ff_image, circles = circles);
    imwrite(v_mean_RGB, fullfile(ToolBox.path_png, sprintf("%s_%s", ToolBox.folder_name, 'v_map.png')))

    % Gif frame generation
    if exportVideos

        parfor frameIdx = 1:numFrames
            v_video_RGB(:, :, :, frameIdx) = flowMapImg(v_rescaled(:, :, frameIdx), {maskArtery, maskVein, maskAV}, {cmapArtery, cmapVein, cmapAV}, background = M0_ff(:, :, frameIdx), circles = circles);
        end

    end

else

    % velocity Colorbars of the gif
    velocityColorbar(cmapArtery, v_min, v_max, 'v_Artery');

    % Average Image
    v_mean_RGB = flowMapImg(v_mean_rescaled, {maskArtery}, {cmapArtery}, background = M0_ff_image, circles = circles);
    imwrite(v_mean_RGB, fullfile(ToolBox.path_png, sprintf("%s_%s", ToolBox.folder_name, 'v_mean.png')))

    % Gif frame generation
    if exportVideos

        parfor frameIdx = 1:numFrames
            v_video_RGB(:, :, :, frameIdx) = flowMapImg(v_rescaled(:, :, frameIdx), {maskArtery}, {cmapArtery}, background = M0_ff(:, :, frameIdx), circles = circles);
        end

    end

end

if exportVideos
    writeGifOnDisc(v_video_RGB, "flowMap");
end

end
