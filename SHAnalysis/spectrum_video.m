function spectrum_video(SH, f_RMS_video, maskArtery, maskBkg)

ToolBox = getGlobalToolBox;
[numX, numY, ~, ~] = size(SH);
f1 = ToolBox.f1;
f2 = ToolBox.f2;
fs = ToolBox.fs;

f_signal = squeeze(sum(f_RMS_video .* maskArtery, [1 2]) / nnz(maskArtery));
f_signal_bkg = squeeze(sum(f_RMS_video .* maskBkg, [1 2]) / nnz(maskBkg));

maskArtery = logical(imresize(maskArtery, [numX, numY]));
maskBkg = logical(imresize(maskBkg, [numX, numY])) & ~maskArtery;

% make video

spectrum_video = spectrum_ploting(SH, f_signal, f_signal_bkg, maskArtery, maskBkg, fs, f1, f2);

writeVideoOnDisc(spectrum_video, fullfile(ToolBox.path_avi, strcat(ToolBox.main_foldername, '_spectrum_video')));
writeGifOnDisc(mat2gray(spectrum_video), "spectrum_video")
end
