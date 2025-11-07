function [maskArtery, maskVein, R_ArterialSignal, diasys_diff] = createMasksSegmentationDeterministic(video, maskVessel, maskArteryTmp)
% createMasksSegmentationDeterministic - Create artery and vein masks using correlation-based segmentation
%   [maskArtery, maskVein, R_ArterialSignal, arterialSignal, quantizedImage, level, color] = createMasksSegmentationDeterministic(video, maskVessel, maskArtery);
%   Inputs:
%   - video: 3D matrix of the video data (height x width x time)
%   - maskVessel: Binary mask for vessels
%   - maskArtery: Binary mask for arteries
% Outputs:
%   - maskArtery: Updated binary mask for arteries
%   - maskVein: Updated binary mask for veins
%   - scoreMaskArtery: Score map for artery mask
%   - scoreMaskVein: Score map for vein mask

ToolBox = getGlobalToolBox;
params = ToolBox.getParams;
mask_params = params.json.Mask;
saveFigures = params.saveFigures;

vesselParams.threshold = mask_params.VascularThreshold;
vesselParams.classes = mask_params.VascularClasses;

M0_ff_img = squeeze(mean(video, 3));

% 1) 0) Compute correlation-based segmentation
[maskArtery_corr, maskVein_corr, R_ArterialSignal, arterialSignal, quantizedImage_corr, level, color] = ...
    correlationSegmentation(video, maskArteryTmp, maskVessel, vesselParams);

if saveFigures
    % 1) 1) Save correlation map
    saveMaskImage(R_ArterialSignal, 'all_21_Correlation.png', isStep = true)
    RGBcorr = labDuoImage(rescale(M0_ff_img), R_ArterialSignal);
    saveMaskImage(RGBcorr, 'all_21_Correlation_rgb.png', isStep = true)

    % 1) 2) Save histogram of correlation values with threshold
    graphThreshHistogram(R_ArterialSignal, level, maskVessel, color, 'all_21_correlation')

    % 1) 3) Save quantized image
    [quantizedImageRGB] = quantizeImageToRGB(quantizedImage_corr, vesselParams.classes);
    saveMaskImage(quantizedImageRGB, 'all_21_quantizedImage_correlation.png', isStep = true, cmap = color);

    t = ToolBox.Cache.t;
    graphSignal('all_21_arterialSignal_correlation', t, squeeze(arterialSignal), '-', color(end, :), ...
        Title = 'Arterial Signal', xlabel = 'Time(s)', ylabel = 'Power Doppler (a.u.)');
end

% 2) 0) Process systolic signal to create artery mask
[M0_Systole_img, M0_Diastole_img] = compute_diasys(video, maskArteryTmp, 'mask');
saveMaskImage(rescale(M0_Systole_img), 'artery_21_systole_img.png', isStep = true)
saveMaskImage(rescale(M0_Diastole_img), 'vein_21_diastole_img.png', isStep = true)

% 2) 1) Compute diastolic-systolic difference image
diasys_diff = M0_Systole_img - M0_Diastole_img;
saveMaskImage(rescale(diasys_diff), 'all_21_diastole_diff.png', isStep = true)
RGBdiasys = labDuoImage(rescale(M0_ff_img), diasys_diff);
saveMaskImage(RGBdiasys, 'vessel_21_diasys_diff.png', isStep = true);

% 2) 2) Create artery and vein masks from diasys signal
[maskArtery_diasys, maskVein_diasys, quantizedImage_diasys] = processDiaSysSignal(diasys_diff, maskVessel, vesselParams, color, 'all_21_diasys');

% 2) 3) Save quantized image
if saveFigures
    [quantizedImageRGB] = quantizeImageToRGB(quantizedImage_diasys, vesselParams.classes);
    saveMaskImage(quantizedImageRGB, 'all_21_quantizedImage_diasys.png', isStep = true, cmap = color);
end

maskArtery = maskArtery_corr | maskArtery_diasys;
maskVein = maskVein_corr | maskVein_diasys;

end
