function [maskArtery, maskVein, maskNeighbors] = createMasks(M0_ff_video, xy_barycenter)
% createMasks - Creates masks for arteries, veins, and neighbors from a video of retinal images.
% Inputs:
%   M0_ff_video: 3D matrix of the video data (height x width x time)
%   xy_barycenter: 2-element vector specifying the barycenter coordinates [x, y]
% Outputs:
%   maskArtery: binary mask for arteries
%   maskVein: binary mask for veins
%   maskNeighbors: binary mask for neighboring pixels

ToolBox = getGlobalToolBox;
params = ToolBox.getParams;
path = ToolBox.path_main;

folder_steps = fullfile('mask', 'steps');

% 0) Initialisation

[numX, numY, numFrames] = size(M0_ff_video);

diaphragmRadius = params.json.Mask.DiaphragmRadius;
cropChoroidRadius = params.json.Mask.CropChoroidRadius;
x_c = xy_barycenter(1);
y_c = xy_barycenter(2);

% Parameters for arteries
vesselParams.threshold = params.json.Mask.VascularThreshold;
vesselParams.classes = params.json.Mask.VascularClasses;

diasysAnalysis = params.json.Mask.DiaSysAnalysis;

% Parameters for arteries
arteryParams.threshold = params.json.Mask.ArterialThreshold;
arteryParams.classes = params.json.Mask.ArterialClasses;

% Parameters for veins
veinParams.threshold = params.json.Mask.VenousThreshold;
veinParams.classes = params.json.Mask.VenousClasses;

minPixelSize = params.json.Mask.MinPixelSize;
forceVesselWidth = params.json.Mask.ForceVesselWidth;

bgWidth = params.json.PulseAnalysis.LocalBackgroundWidth;
r1 = params.json.SizeOfField.SmallRadiusRatio;
r2 = params.json.SizeOfField.BigRadiusRatio;

cArtery = [0, 0, 0; 255, 22, 18] / 255;
cVein = [0, 0, 0; 18, 23, 255] / 255;

cmapArtery = cmapLAB(256, [0 0 0], 0, [1 0 0], 1/3, [1 1 0], 2/3, [1 1 1], 1);
cmapVein = cmapLAB(256, [0 0 0], 0, [0 0 1], 1/3, [0 1 1], 2/3, [1 1 1], 1);
cmapAV = cmapLAB(256, [0 0 0], 0, [1 0 1], 1/3, [1 1 1], 1);

% 1) First Masks and Correlation

M0_ff_img = squeeze(mean(M0_ff_video, 3));
saveImage(M0_ff_img, 'all_10_M0.png', isStep = true)

maskDiaphragm = diskMask(numX, numY, diaphragmRadius);
saveImage(rescale(M0_ff_img) + maskDiaphragm .* 0.5, 'all_11_maskDiaphragm.png', isStep = true)

M0_video = M0_ff_video;
A = ones(1, 1, numFrames);
B = A .* maskDiaphragm;
M0_video(~B) = NaN;
M0_img = squeeze(mean(M0_video, 3, 'omitnan'));

% 1) 1) Compute vesselness response
% Frangi and Gabor Vesselness

% Frangi Vesselness
[maskVesselnessFrangi, M0_Frangi] = frangiVesselness(M0_img, ...
    'range', params.json.Mask.VesselnessFrangiRange, ...
    'step', params.json.Mask.VesselnessFrangiStep);

saveImage(maskVesselnessFrangi, 'all_12_frangi_mask.png', isStep = true)
saveImage(M0_Frangi, 'all_12_frangi_img.png', isStep = true)

% Gabor Vesselness
[maskVesselnessGabor, M0_Gabor] = gaborVesselness(M0_ff_img, ...
    'range', params.json.Mask.VesselnessGaborRange, ...
    'step', params.json.Mask.VesselnessGaborStep);

saveImage(maskVesselnessGabor & maskDiaphragm, 'all_12_gabor_mask.png', isStep = true)
saveImage(M0_Gabor, 'all_12_gabor_img.png', isStep = true)

if params.json.Mask.VesselnessHolonet

    try
        maskVesselness = getHolonetprediction(M0_ff_img);
        saveImage(maskVesselness, 'all_14_maskHoloNet.png', isStep = true)
    catch
        warning("The Holonet ONNX-model couldn't be found.")
        maskVesselness = (maskVesselnessFrangi | maskVesselnessGabor) & maskDiaphragm;
        saveImage(maskVesselness, 'all_14_maskDeterministic.png', isStep = true)
    end

else
    maskVesselness = (maskVesselnessFrangi | maskVesselnessGabor) & maskDiaphragm;
    saveImage(maskVesselness, 'all_14_maskDeterministic.png', isStep = true)
end

% 1) 2) Compute the barycenters and the circle mask

maskCircle = diskMask(numX, numY, cropChoroidRadius, 'center', [x_c / numX, y_c / numY]);
maskVesselnessClean = maskVesselness & bwareafilt(maskVesselness | maskCircle, 1, 8);
saveImage(maskVesselnessClean + maskCircle * 0.5, 'all_15_VesselMask_clear.png', isStep = true)

%  1) 3) Compute first correlation

cVascular = [0 0 0];

if ~isempty(vesselParams.threshold)

    if abs(vesselParams.threshold) <= 1
        vesselParams.classes = [-1; 1];
    end

end

[maskArtery, maskVein, R_VascularSignal, vascularSignal, quantizedImage, level, color] = correlationSegmentation(M0_video, maskVesselnessClean, vesselParams);

% 1) 4) Save all images

t = linspace(0, numFrames * ToolBox.stride / ToolBox.fs / 1000, numFrames);
graphSignal('all_15_vascularSignal', folder_steps, t, squeeze(vascularSignal), '-', cVascular, ...
    Title = 'Vascular Signal', xlabel = 'Time(s)', ylabel = 'Power Doppler (a.u.)');

saveImage(R_VascularSignal, 'all_15_Correlation.png', isStep = true)
RGBcorr = labDuoImage(M0_Gabor, R_VascularSignal);
saveImage(RGBcorr, 'all_15_Correlation_rgb.png', isStep = true)

[quantizedImageRGB] = quantizeImageToRGB(quantizedImage, vesselParams.classes);
saveImage(quantizedImageRGB, 'all_16_quantizedImage.png', isStep = true, cmap = color);
graphThreshHistogram(R_VascularSignal, level, maskVesselnessClean, color, 'all_16')

saveImage(maskArtery, 'artery_17_FirstMask.png', isStep = true, cmap = cArtery)
saveImage(maskVein, 'vein_17_FirstMask.png', isStep = true, cmap = cVein)

% 1) 5) Remove small blobs
maskArteryClean = bwareaopen(maskArtery, minPixelSize);
maskVeinClean = bwareaopen(maskVein, minPixelSize);

saveImage(maskArteryClean, 'artery_18_FirstMaskClean.png', isStep = true, cmap = cArtery)
saveImage(maskVeinClean, 'vein_18_FirstMaskClean.png', isStep = true, cmap = cVein)

M0_Artery = setcmap(M0_ff_img, maskArtery, cmapArtery);
M0_Vein = setcmap(M0_ff_img, maskVein, cmapVein);
M0_AV = setcmap(M0_ff_img, maskArtery & maskVein, cmapAV);

M0_RGB = (M0_Artery + M0_Vein) .* ~(maskArtery & maskVein) + M0_AV + rescale(M0_ff_img) .* ~(maskArtery | maskVein);
saveImage(M0_RGB, 'all_19_RGB.png', isStep = true)

% 2)  Improvements of the first mask

if params.json.Mask.ImproveMask

    % 2) 0) Computation of the M0 in Diastole and in Systole

    [M0_Systole_img, M0_Diastole_img, M0_Systole_video] = compute_diasys(M0_video, maskArtery, 'mask');
    [M0_Sys_img, M0_Dia_img, ~] = compute_diasys(M0_ff_video, maskArtery, 'mask');
    saveImage(rescale(M0_Systole_img), 'artery_20_systole_img.png', isStep = true)
    saveImage(rescale(M0_Diastole_img), 'vein_20_diastole_img.png', isStep = true)

    % 2) 1) New Vesselness Mask

    Systole_Frangi = frangiVesselness(M0_Systole_img, ...
        'range', params.json.Mask.VesselnessFrangiRange, ...
        'step', params.json.Mask.VesselnessFrangiStep);
    Diastole_Frangi = frangiVesselness(M0_Diastole_img, ...
        'range', params.json.Mask.VesselnessFrangiRange, ...
        'step', params.json.Mask.VesselnessFrangiStep);

    saveImage(Systole_Frangi, 'artery_20_frangi_mask.png', isStep = true)
    saveImage(Diastole_Frangi, 'vein_20_frangi_mask.png', isStep = true)

    Systole_Gabor = gaborVesselness(M0_Sys_img, ...
        'range', params.json.Mask.VesselnessGaborRange, ...
        'step', params.json.Mask.VesselnessGaborStep);
    Diastole_Gabor = gaborVesselness(M0_Dia_img, ...
        'range', params.json.Mask.VesselnessGaborRange, ...
        'step', params.json.Mask.VesselnessGaborStep);

    saveImage(Systole_Gabor & maskDiaphragm, 'artery_20_gabor_mask.png', isStep = true)
    saveImage(Diastole_Gabor & maskDiaphragm, 'vein_20_gabor_mask.png', isStep = true)

    if params.json.Mask.VesselnessHolonet
        maskVesselness = maskVesselness & maskDiaphragm;
    else
        maskVesselness = (Systole_Frangi | Diastole_Frangi | Systole_Gabor | Diastole_Gabor) & maskDiaphragm;
    end

    maskVesselnessClean = maskVesselness & bwareafilt(maskVesselness | maskCircle, 1, 8);
    saveImage(maskVesselnessClean + maskCircle * 0.5, 'all_20_VesselMask_clear.png', isStep = true)

    % 2) 2) Diastole-Systole Image

    diasysArtery = M0_Systole_img - M0_Diastole_img;
    mDiasys = mean(diasysArtery, 'all', 'omitnan');
    diasysVein = mDiasys - diasysArtery;
    saveImage(diasysArtery, 'artery_21_diasys_img.png', isStep = true)
    saveImage(diasysVein, 'vein_21_diasys_img.png', isStep = true)

    RGBdiasys = labDuoImage(rescale(M0_Gabor), diasysArtery);
    saveImage(RGBdiasys, 'vessel_40_diasys_rgb.png', isStep = true)
    saveImage(RGBdiasys, 'DiaSysRGB.png')

    if ~isempty(arteryParams.threshold)

        if abs(arteryParams.threshold) <= 1
            arteryParams.classes = [-1; 1];
        end

    end

    if ~isempty(veinParams.threshold)

        if abs(veinParams.threshold) <= 1
            veinParams.classes = [-1; 1];
        end

    end

    if diasysAnalysis % Systole/Diastole Analysis

        % Artery
        [maskArtery, ~, quantizedImage] = processDiaSysSignal(diasysArtery, maskVesselnessClean, arteryParams, cArtery, 'artery_23');
        saveImage(maskArtery, 'artery_23_DiaSysMask.png', isStep = true, cmap = cArtery);
        saveImage(quantizedImage, 'artery_23_quantizedImage.png', isStep = true);

        quantizedImageRGB = quantizeImageToRGB(quantizedImage, arteryParams.classes);
        saveImage(quantizedImageRGB, 'artery_23_quantizedImage.png', isStep = true, cmap = cArtery);

        % Vein
        [~, maskVein, quantizedImage] = processDiaSysSignal(diasysVein, maskVesselnessClean, veinParams, cVein, 'vein_23');
        saveImage(maskVein, 'vein_23_DiaSysMask.png', isStep = true, cmap = cVein);
        saveImage(quantizedImage, 'vein_23_quantizedImage.png', isStep = true);

        quantizedImageRGB = quantizeImageToRGB(quantizedImage, veinParams.classes);
        saveImage(quantizedImageRGB, 'vein_23_quantizedImage.png', isStep = true, cmap = cVein);
    else % Second Correlation Analysis

        % Artery
        [maskArtery, ~, R, ~, quantizedImage, level, color] = correlationSegmentation(M0_Systole_video, maskVesselnessClean, arteryParams);
        saveImage(R, 'artery_23_Correlation.png', isStep = true)

        RGBcorr = labDuoImage(M0_Gabor, R);
        saveImage(RGBcorr, 'artery_23_Correlation_rgb.png', isStep = true)

        graphThreshHistogram(R, level, maskVesselnessClean, color, 'artery_23');
        saveImage(quantizedImage, 'artery_23_quantizedImage.png', isStep = true);

        quantizedImageRGB = quantizeImageToRGB(quantizedImage, arteryParams.classes);
        saveImage(quantizedImageRGB, 'artery_23_quantizedImage_rgb.png', isStep = true, cmap = color);

        % Vein
        [~, maskVein, quantizedImage] = processDiaSysSignal(diasysVein, maskVesselnessClean, veinParams, cVein, 'vein_23');
        saveImage(maskVein, 'vein_23_DiaSysMask.png', isStep = true, cmap = cVein);
        saveImage(quantizedImage, 'vein_23_quantizedImage.png', isStep = true);

        quantizedImageRGB = quantizeImageToRGB(quantizedImage, veinParams.classes);
        saveImage(quantizedImageRGB, 'vein_23_quantizedImage_rgb.png', isStep = true);

    end

    % 3) Mask Clearing

    % 3) 0) Morphological Operations
    [mask_dilated, mask_closed, mask_opened, mask_widened] = clearMasks(maskArtery, ...
        'min_area', params.json.Mask.MinPixelSize, ...
        'imclose_radius', params.json.Mask.ImcloseRadius, ...
        'min_width', params.json.Mask.MinimumVesselWidth, ...
        'imdilate_size', params.json.Mask.FinalDilation);

    saveImage(mask_dilated, 'artery_30_ClearedMask.png', isStep = true, cmap = cArtery);
    saveImage(mask_opened, 'artery_30_ClearedMask_opened.png', isStep = true, cmap = cArtery);
    saveImage(mask_closed, 'artery_30_ClearedMask_closed.png', isStep = true, cmap = cArtery);
    saveImage(mask_widened, 'artery_30_ClearedMask_widened.png', isStep = true, cmap = cArtery);

    maskArtery = mask_dilated;

    [mask_dilated, mask_closed, mask_opened, mask_widened] = clearMasks(maskVein, ...
        'min_area', params.json.Mask.MinPixelSize, ...
        'imclose_radius', params.json.Mask.ImcloseRadius, ...
        'min_width', params.json.Mask.MinimumVesselWidth, ...
        'imdilate_size', params.json.Mask.FinalDilation);

    saveImage(mask_dilated, 'vein_30_ClearedMask.png', isStep = true, cmap = cVein);
    saveImage(mask_opened, 'vein_30_ClearedMask_opened.png', isStep = true, cmap = cVein);
    saveImage(mask_closed, 'vein_30_ClearedMask_closed.png', isStep = true, cmap = cVein);
    saveImage(mask_widened, 'vein_30_ClearedMask_widened.png', isStep = true, cmap = cVein);

    maskVein = mask_dilated;

    % 3) 1) Final Blob removal
    maskVessel = maskArtery | maskVein;

    maskArtery = maskArtery & bwareafilt(maskVessel | maskCircle, 1, 8);
    saveImage(maskArtery + maskCircle * 0.5, 'artery_31_VesselMask_clear.png', isStep = true)

    maskVein = maskVein & bwareafilt(maskVessel | maskCircle, 1, 8);
    saveImage(maskVein + maskCircle * 0.5, 'vein_31_VesselMask_clear.png', isStep = true)

    maskArtery_no_import = maskArtery;
    maskVein_no_import = maskVein;

    % 3) 2) Force Create Masks in case they exist
    if (params.json.Mask.ForcedMasks == -1 || params.json.Mask.ForcedMasks == 1)

        if isfile(fullfile(path, 'mask', 'forceMaskArtery.png'))
            maskArtery = mat2gray(mean(imread(fullfile(path, 'mask', 'forceMaskArtery.png')), 3)) > 0;

            if size(maskArtery, 1) ~= maskCircle
                maskArtery = imresize(maskArtery, [numX, numY], "nearest");
            end

        elseif params.json.Mask.ForcedMasks == 1
            error("Cannot force Artery Mask because none given in the mask folder. Please create a forceMaskArtery.png file in the mask folder. (SET ForcedMasks to -1 to skip use the auto mask)");
        end

        if isfile(fullfile(path, 'mask', 'forceMaskVein.png'))
            maskVein = mat2gray(mean(imread(fullfile(path, 'mask', 'forceMaskVein.png')), 3)) > 0;

            if size(maskVein, 1) ~= maskCircle
                maskVein = imresize(maskVein, [numX, numY], "nearest");
            end

        elseif params.json.Mask.ForcedMasks == 1
            error("Cannot force Vein Mask because none given in the mask folder. Please create a forceMaskVein.png file in the mask folder. (SET ForcedMasks to -1 to skip use the auto mask)");
        end

    end

    % 3) 3) Segmentation Scores Calculation

    segmentationScores(maskArtery, maskVein);

    % 3) 4) Segmention force width

    if forceVesselWidth > 0
        dilationSE = strel('disk', forceVesselWidth);
        maskArtery = imdilate(bwskel(maskArtery), dilationSE);
        maskVein = imdilate(bwskel(maskVein), dilationSE);

        maskArtery = imdilate(bwskel(maskArtery), dilationSE);
        maskVein = imdilate(bwskel(maskVein), dilationSE);
    end

else

    [mask_dilated, mask_closed, mask_opened, mask_widened] = clearMasks(maskArtery, ...
        'min_area', params.json.Mask.MinPixelSize, ...
        'imclose_radius', params.json.Mask.ImcloseRadius, ...
        'min_width', params.json.Mask.MinimumVesselWidth, ...
        'imdilate_size', params.json.Mask.FinalDilation);

    saveImage(mask_dilated, 'artery_30_ClearedMask.png', isStep = true, cmap = cArtery);
    saveImage(mask_opened, 'artery_30_ClearedMask_opened.png', isStep = true, cmap = cArtery);
    saveImage(mask_closed, 'artery_30_ClearedMask_closed.png', isStep = true, cmap = cArtery);
    saveImage(mask_widened, 'artery_30_ClearedMask_widened.png', isStep = true, cmap = cArtery);

    maskArtery = mask_dilated;

    [mask_dilated, mask_closed, mask_opened, mask_widened] = clearMasks(maskVein, ...
        'min_area', params.json.Mask.MinPixelSize, ...
        'imclose_radius', params.json.Mask.ImcloseRadius, ...
        'min_width', params.json.Mask.MinimumVesselWidth, ...
        'imdilate_size', params.json.Mask.FinalDilation);

    saveImage(mask_dilated, 'vein_30_ClearedMask.png', isStep = true, cmap = cVein);
    saveImage(mask_opened, 'vein_30_ClearedMask_opened.png', isStep = true, cmap = cVein);
    saveImage(mask_closed, 'vein_30_ClearedMask_closed.png', isStep = true, cmap = cVein);
    saveImage(mask_widened, 'vein_30_ClearedMask_widened.png', isStep = true, cmap = cVein);

    maskVein = mask_dilated;

    maskVessel = maskArtery | maskVein;

    maskArtery = maskArtery & bwareafilt(maskVessel | maskCircle, 1, 8);
    saveImage(maskArtery + maskCircle * 0.5, 'artery_31_VesselMask_clear.png', isStep = true)

    maskVein = maskVein & bwareafilt(maskVessel | maskCircle, 1, 8);
    saveImage(maskVein + maskCircle * 0.5, 'vein_31_VesselMask_clear.png', isStep = true)

    maskArtery_no_import = maskArtery;
    maskVein_no_import = maskVein;

    % 3) 2) Force Create Masks in case they exist
    if (params.json.Mask.ForcedMasks == -1 || params.json.Mask.ForcedMasks == 1)

        if isfile(fullfile(path, 'mask', 'forceMaskArtery.png'))
            maskArtery = mat2gray(mean(imread(fullfile(path, 'mask', 'forceMaskArtery.png')), 3)) > 0;

            if size(maskArtery, 1) ~= maskCircle
                maskArtery = imresize(maskArtery, [numX, numY], "nearest");
            end

        elseif params.json.Mask.ForcedMasks == 1
            error("Cannot force Artery Mask because none given in the mask folder. Please create a forceMaskArtery.png file in the mask folder. (SET ForcedMasks to -1 to skip use the auto mask)");
        end

        if isfile(fullfile(path, 'mask', 'forceMaskVein.png'))
            maskVein = mat2gray(mean(imread(fullfile(path, 'mask', 'forceMaskVein.png')), 3)) > 0;

            if size(maskVein, 1) ~= maskCircle
                maskVein = imresize(maskVein, [numX, numY], "nearest");
            end

        elseif params.json.Mask.ForcedMasks == 1
            error("Cannot force Vein Mask because none given in the mask folder. Please create a forceMaskVein.png file in the mask folder. (SET ForcedMasks to -1 to skip use the auto mask)");
        end

    end

end

% 3) 5) Create Vessel and Background Mask
maskArtery = maskArtery & maskDiaphragm;
maskVein = maskVein & maskDiaphragm;
maskVessel = maskArtery | maskVein;
maskBackground = not(maskVessel);

% 4) FINAL FIGURES

% 4) 1) RGB Figures
if isfile(fullfile(path, 'mask', 'forceMaskArtery.png')) || isfile(fullfile(path, 'mask', 'forceMaskVein.png'))

    maskVesselNoImport = maskArtery_no_import | maskVein_no_import;

    maskArtery_no_import = maskArtery_no_import & bwareafilt(maskVesselNoImport | maskCircle, 1, 8);
    saveImage(maskArtery_no_import + maskCircle * 0.5, 'artery_31_VesselMask_clear.png', isStep = true)

    maskVein_no_import = maskVein_no_import & bwareafilt(maskVesselNoImport | maskCircle, 1, 8);
    saveImage(maskVein_no_import + maskCircle * 0.5, 'vein_31_VesselMask_clear.png', isStep = true)

    M0_Artery = setcmap(M0_ff_img, maskArtery_no_import, cmapArtery);
    M0_Vein = setcmap(M0_ff_img, maskVein_no_import, cmapVein);
    M0_AV = setcmap(M0_ff_img, maskArtery_no_import & maskVein_no_import, cmapAV);
    M0_RGB = (M0_Artery + M0_Vein) .* ~(maskArtery_no_import & maskVein_no_import) + M0_AV + rescale(M0_ff_img) .* ~(maskArtery_no_import | maskVein_no_import);
    saveImage(M0_RGB, 'vessel_40_RGB_no_import.png', isStep = true)
end

M0_Artery = setcmap(M0_ff_img, maskArtery, cmapArtery);
M0_Vein = setcmap(M0_ff_img, maskVein, cmapVein);
M0_AV = setcmap(M0_ff_img, maskArtery & maskVein, cmapAV);

M0_RGB = (M0_Artery + M0_Vein) .* ~(maskArtery & maskVein) + M0_AV + rescale(M0_ff_img) .* ~(maskArtery | maskVein);
saveImage(M0_RGB, 'vessel_40_RGB.png', isStep = true)
saveImage(M0_RGB, 'RGB_img.png')

% 4) 2) Neighbours Mask

if params.json.Mask.AllNonVesselsAsBackground
    maskNeighbors = (maskBackground & ~maskVesselness) & maskDiaphragm;
else
    if params.veins_analysis
        maskNeighbors = imdilate(maskArtery | maskVein, strel('disk', bgWidth)) & ~(maskArtery | maskVein);
    else
        maskNeighbors = imdilate(maskArtery, strel('disk', bgWidth)) & ~(maskArtery); % & ~(maskArtery | maskVein); possible 
    end
end

cmapNeighbors = cmapLAB(256, [0 1 0], 0, [1 1 1], 1);

M0_Neighbors = setcmap(M0_ff_img, maskNeighbors, cmapNeighbors);


if params.veins_analysis
    neighborsMaskSeg = (M0_Artery + M0_Vein) .* ~(maskArtery & maskVein) + ...
        M0_AV + M0_Neighbors + ...
        rescale(M0_ff_img) .* ~(maskArtery | maskVein | maskNeighbors);
else
    neighborsMaskSeg = (M0_Artery) .* ~(maskArtery) + ...
        M0_Artery + M0_Neighbors + ...
        rescale(M0_ff_img) .* ~(maskArtery | maskNeighbors);
end

saveImage(neighborsMaskSeg, 'neighbors_img.png')

% 4) 4) Save all images

saveImage(maskArtery, 'maskArtery.png')
saveImage(maskVein, 'maskVein.png')
saveImage(maskVessel, 'maskVessel.png')
saveImage(maskNeighbors, 'maskNeighbors.png')
saveImage(maskBackground, 'maskBackground.png')
saveImage(bwskel(maskArtery), 'skeletonArtery.png')
saveImage(bwskel(maskVein), 'skeletonVein.png')

% 4) 5) Mask Section & Force Barycenter

xy_barycenter = [x_c, y_c];
createMaskSection(ToolBox, M0_ff_img, r1, r2, xy_barycenter, 'vesselMapArtery', maskArtery, thin = 0.01);
createMaskSection(ToolBox, M0_ff_img, r1, r2, xy_barycenter, 'vesselMap', maskArtery, maskVein, thin = 0.01);

% 4) 6) Arteries tree
% try
%     getLongestArteryBranch(maskArtery, xy_barycenter, 'Artery');
%     getLongestArteryBranch(maskVein, xy_barycenter, 'Vein');
% catch ME
%
%     for i = 1:length(ME.stack)
%         disp("Error in getLongestArteryBranch: " + ME.stack(i).name + " at line " + ME.stack(i).line)
%     end
%
% end

close all
end
