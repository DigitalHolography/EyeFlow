function [v_video_RGB, v_mean_RGB] = bloodFlowVelocity(v_video, maskArtery, maskVein, M0_ff_video, xy_barycenter)

close all
ToolBox = getGlobalToolBox;
params = ToolBox.getParams;
veinsAnalysis = params.veins_analysis;
exportVideos = params.exportVideos;

% Rescale once
M0_ff_video = rescale(M0_ff_video);
[numX, numY, numFrames] = size(v_video);
x_c = xy_barycenter(1) / numX;
y_c = xy_barycenter(2) / numY;
r1 = params.json.SizeOfField.SmallRadiusRatio;
r2 = params.json.SizeOfField.BigRadiusRatio;
maskSection = diskMask(numX, numY, r1, r2, center = [x_c, y_c]);

% Precompute masks
maskAV = maskArtery & maskVein;
maskArterySection = maskArtery & maskSection & ~maskAV;
maskVeinSection = maskVein & maskSection & ~maskAV;

% 1) VELOCITY VIDEO
tVelocityVideo = tic;

[v_video_RGB, v_mean_RGB] = flowMap(v_video, maskSection, maskArtery, maskVein, M0_ff_video, xy_barycenter, ToolBox);

fprintf("- Velocity Map Timing : %ds\n", round(toc(tVelocityVideo)))

% 2) HISTOGRAM

histoVideoArtery = VelocityHistogram(v_video, maskArterySection, 'Artery');

if veinsAnalysis
    histoVideoVein = VelocityHistogram(v_video, maskVeinSection, 'Vein');
end

% 3) COMBINED

[numX_fig, numY_fig, ~, ~] = size(histoVideoArtery);

if veinsAnalysis
    v_mean_RGB4Gif = rescale(imresize3(v_mean_RGB, [numX_fig*2 numX_fig*2 3]));
else
    v_mean_RGB4Gif = rescale(imresize3(v_mean_RGB, [numY_fig numY_fig 3]));
end

if exportVideos

    if veinsAnalysis
        imwrite(cat(2, v_mean_RGB4Gif, cat(1, mat2gray(histoVideoArtery(:, :, :, end)), mat2gray(histoVideoVein(:, :, :, end)))), fullfile(ToolBox.path_png, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, 'AVGflowVideoCombined.png')))
    else
        imwrite(cat(1, v_mean_RGB4Gif, mat2gray(histoVideoArtery(:, :, :, end))), fullfile(ToolBox.path_png, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, 'AVGflowVideoCombined.png')))
    end

    if veinsAnalysis
        v_video_RGB4Gif(:, :, 1, :) = mat2gray(imresize3(squeeze(v_video_RGB(:, :, 1, :)), [numX_fig*2 numX_fig*2 numFrames]));
        v_video_RGB4Gif(:, :, 2, :) = mat2gray(imresize3(squeeze(v_video_RGB(:, :, 2, :)), [numX_fig*2 numX_fig*2 numFrames]));
        v_video_RGB4Gif(:, :, 3, :) = mat2gray(imresize3(squeeze(v_video_RGB(:, :, 3, :)), [numX_fig*2 numX_fig*2 numFrames]));
        combinedGifs = cat(2, v_video_RGB4Gif, cat(1, mat2gray(histoVideoArtery), mat2gray(histoVideoVein)));
    else
        v_video_RGB4Gif(:, :, 1, :) = mat2gray(imresize3(squeeze(v_video_RGB(:, :, 1, :)), [numY_fig numY_fig numFrames]));
        v_video_RGB4Gif(:, :, 2, :) = mat2gray(imresize3(squeeze(v_video_RGB(:, :, 2, :)), [numY_fig numY_fig numFrames]));
        v_video_RGB4Gif(:, :, 3, :) = mat2gray(imresize3(squeeze(v_video_RGB(:, :, 3, :)), [numY_fig numY_fig numFrames]));
        combinedGifs = cat(1, v_video_RGB4Gif, mat2gray(histoVideoArtery));
    end

    writeGifOnDisc(mat2gray(combinedGifs), "velocityHistogramCombined", 0.04);

else

    if veinsAnalysis
        imwrite(cat(2, v_mean_RGB4Gif, cat(1, mat2gray(histoVideoArtery), mat2gray(histoVideoVein))), fullfile(ToolBox.path_png, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, 'AVGflowVideoCombined.png')))
    else
        imwrite(cat(1, v_mean_RGB4Gif, mat2gray(histoVideoArtery)), fullfile(ToolBox.path_png, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, 'AVGflowVideoCombined.png')))
    end

end

close all

end
