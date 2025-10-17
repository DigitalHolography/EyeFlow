function [] = spectrum_analysis(SH_cube, M0_ff)

ToolBox = getGlobalToolBox;

fs = ToolBox.fs / 2;
f1 = ToolBox.f1;
f2 = (ToolBox.f1 + ToolBox.f2) / 2;
f3 = ToolBox.f2;

[~, ~, numFreq, numFrames] = size(SH_cube);
[numX, numY, ~] = size(M0_ff);

% integration intervals
low_n1 = round(f1 * numFreq / fs) + 1;
low_n2 = floor(f2 * numFreq / fs);
high_n1 = ceil(f2 * numFreq / fs);
high_n2 = round(f3 * numFreq / fs);

I = rescale(mean(M0_ff, 3));
MeanFreqLow = sum(SH_cube(:, :, low_n1:low_n2, :), [3 4]) / (low_n2 - low_n1) / numFrames;
MeanFreqHigh = sum(SH_cube(:, :, high_n1:high_n2, :), [3 4]) / (high_n2 - high_n1) / numFrames;

df = imresize(MeanFreqHigh - MeanFreqLow, [numX numY]);
m = sum(df .* fftshift(diskMask(numX, numY, 0.1)), [1 2]) ./ nnz(fftshift(diskMask(numX, numY, 0.1)));
I_lab = labDuoImage(I, (df - m));

imwrite(I_lab, fullfile(ToolBox.path_png, 'spectralAnalysis', sprintf("%s_ColorImg.png", ToolBox.folder_name)))

SH_ColorVideoRGB = zeros(numX, numY, 3, numFrames);

parfor frameIdx = 1:numFrames

    % integration
    freq_low = squeeze(sum(abs(SH_cube(:, :, low_n1:low_n2, frameIdx)), 3)) / (low_n2 - low_n1);
    freq_high = squeeze(sum(abs(SH_cube(:, :, high_n1:high_n2, frameIdx)), 3)) / (high_n2 - high_n1);
    df = imresize(freq_high - freq_low, [numX numY]);
    m = sum(df .* fftshift(diskMask(numX, numY, 0.1)), [1 2]) ./ nnz(fftshift(diskMask(numX, numY, 0.1)));
    I_lab = labDuoImage(I, (df - m));
    SH_ColorVideoRGB(:, :, :, frameIdx) = mat2gray(I_lab);

end

% save video
writeVideoOnDisc(SH_ColorVideoRGB, fullfile(ToolBox.path_avi, strcat(ToolBox.folder_name, '_SH_ColorVideo')));
writeGifOnDisc(SH_ColorVideoRGB, 'ColorVideo');

SH_ColorSpectrumRGB = zeros(numX, numY, 3, numFreq - 2);

parfor freqIdx = low_n1 + 1:high_n2 - 1

    % integration
    freq_low = squeeze(sum(abs(SH_cube(:, :, low_n1:freqIdx, :)), [3 4])) / (freqIdx - low_n1) / numFrames;
    freq_high = squeeze(sum(abs(SH_cube(:, :, freqIdx:high_n2, :)), [3 4])) / (high_n2 - freqIdx) / numFrames;
    df = imresize(rescale(freq_high) - rescale(freq_low), [numX numY]);
    m = sum(df .* fftshift(diskMask(numX, numY, 0.1)), [1 2]) ./ nnz(fftshift(diskMask(numX, numY, 0.1)));
    I_lab = labDuoImage(I, (df - m), 180);
    SH_ColorSpectrumRGB(:, :, :, freqIdx) = mat2gray(I_lab);

end

% save video
writeVideoOnDisc(SH_ColorSpectrumRGB, fullfile(ToolBox.path_avi, strcat(ToolBox.folder_name, '_SH_ColorSpectra')));
writeGifOnDisc(SH_ColorSpectrumRGB(:, :, :, low_n1 + 1:high_n2 - 1), 'ColorSpectra', 0.1, high_n2 - low_n1 - 1);

end
