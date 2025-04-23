function spectrum_video(SH, f_RMS_video, maskArtery, maskBkg)

ToolBox = getGlobalToolBox;
[numX, numY, ~, numFrames] = size(SH);
f1 = ToolBox.f1;
f2 = ToolBox.f2;
fs = ToolBox.fs;

f_signal = squeeze(sum(f_RMS_video .* maskArtery, [1 2]) / nnz(maskArtery));
f_signal_bkg = squeeze(sum(f_RMS_video .* maskBkg, [1 2]) / nnz(maskBkg));

maskArteryResized = logical(imresize(maskArtery, [numX, numY]));
maskBkgResized = logical(imresize(maskBkg, [numX, numY])) & ~maskArteryResized;

spectrum_video = zeros(420, 560, 3, numFrames);
% make video
parfor frameIdx = 1:numFrames
    fi = figure("Visible", "off", "Color", 'w');
    spectrum_ploting(SH(:, :, :, frameIdx), f_signal(frameIdx), f_signal_bkg(frameIdx), ...
        maskArteryResized, maskBkgResized, fs, f1, f2);

    frame = getframe(fi);
    spectrum_video(:, :, :, frameIdx) = frame.cdata;
    close(fi)
end

writeVideoOnDisc(spectrum_video, fullfile(ToolBox.path_avi, strcat(ToolBox.main_foldername, '_spectrum_video')));
writeGifOnDisc(mat2gray(spectrum_video), "spectrum_video.gif")
end
