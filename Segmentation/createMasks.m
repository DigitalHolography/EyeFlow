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
cVessels = [18, 23, 255; 255, 22, 18] / 255;

cmapArtery = cmapLAB(256, [0 0 0], 0, [1 0 0], 1/3, [1 1 0], 2/3, [1 1 1], 1);
cmapVein = cmapLAB(256, [0 0 0], 0, [0 0 1], 1/3, [0 1 1], 2/3, [1 1 1], 1);
cmapAV = cmapLAB(256, [0 0 0], 0, [1 0 1], 1/3, [1 1 1], 1);

% 1) First Masks and Correlation

M0_ff_img = squeeze(mean(M0_ff_video, 3));
saveImage(M0_ff_img, ToolBox, 'all_10_M0.png', isStep = true)

if ~isfile(fullfile(ToolBox.path_gif, sprintf("%s_M0.gif", ToolBox.folder_name)))
    writeGifOnDisc(imresize(rescale(M0_ff_video), 0.5), "M0")
end

maskDiaphragm = diskMask(numX, numY, diaphragmRadius);
saveImage(rescale(M0_ff_img) + maskDiaphragm .* 0.5, ToolBox, 'all_11_maskDiaphragm.png', isStep = true)

M0_video = M0_ff_video;
A = ones(1, 1, numFrames);
B = A .* maskDiaphragm;
M0_video(~B) = NaN;
M0_img = squeeze(mean(M0_video, 3, 'omitnan'));

% 1) 1) Compute vesselness response

[maskVesselnessFrangi] = frangiVesselness(M0_img, 'all_12', ToolBox);
[maskVesselnessGabor, M0_Gabor] = gaborVesselness(M0_ff_img, ToolBox, 'all_13');

if params.json.Mask.VesselnessHolonet

    try
        maskVesselness = getHolonetprediction(M0_ff_img);
        saveImage(maskVesselness, ToolBox, 'all_14_maskHoloNet.png', isStep = true)
    catch
        warning("The Holonet ONNX-model couldn't be found.")
        maskVesselness = (maskVesselnessFrangi | maskVesselnessGabor) & maskDiaphragm;
        saveImage(maskVesselness, ToolBox, 'all_14_maskDeterministic.png', isStep = true)
    end

else
    maskVesselness = (maskVesselnessFrangi | maskVesselnessGabor) & maskDiaphragm;
    saveImage(maskVesselness, ToolBox, 'all_14_maskDeterministic.png', isStep = true)
end

% 1) 2) Compute the barycenters and the circle mask

maskCircle = diskMask(numX, numY, cropChoroidRadius, 'center', [x_c / numX, y_c / numY]);
maskVesselnessClean = removeDisconnected(maskVesselness, maskVesselness, maskCircle, 'all_15_VesselMask', ToolBox);

%  1) 3) Compute first correlation

cVascular = [0 0 0];
[maskArtery, maskVein, R_VascularSignal, vascularSignal, quantizedImage, level, color] = correlationSegmentation(M0_video, maskVesselnessClean, vesselParams);

% 1) 4) Save all images

t = linspace(0, numFrames * ToolBox.stride / ToolBox.fs / 1000, numFrames);
graphSignal('all_15_vascularSignal', folder_steps, t, squeeze(vascularSignal), '-', cVascular, ...
    Title = 'Vascular Signal', xlabel = 'Time(s)', ylabel = 'Power Doppler (a.u.)');

saveImage(R_VascularSignal, ToolBox, 'all_15_Correlation.png', isStep = true)
RGBcorr = labDuoImage(M0_ff_img, R_VascularSignal);
saveImage(RGBcorr, ToolBox, 'all_15_Correlation_rgb.png', isStep = true)

if ~isempty(vesselParams.threshold)
    % IF Manual Thresholds have been set between -1 and 1 then they are used
    graphThreshHistogram(R_VascularSignal, vesselParams.threshold, maskVesselnessClean, cVessels, 'all_16')
else
    [quantizedImageRGB] = quantizeImageToRGB(quantizedImage, vesselParams.classes);
    saveImage(quantizedImageRGB, ToolBox, 'all_16_quantizedImage.png', isStep = true, cmap = color);
    graphThreshHistogram(R_VascularSignal, level, maskVesselnessClean, color, 'all_16');

end

saveImage(maskArtery, ToolBox, 'artery_17_FirstMask.png', isStep = true, cmap = cArtery)
saveImage(maskVein, ToolBox, 'vein_17_FirstMask.png', isStep = true, cmap = cVein)

% 1) 5) Remove small blobs
maskArteryClean = bwareaopen(maskArtery, minPixelSize);
maskVeinClean = bwareaopen(maskVein, minPixelSize);

saveImage(maskArteryClean, ToolBox, 'artery_18_FirstMaskClean.png', isStep = true, cmap = cArtery)
saveImage(maskVeinClean, ToolBox, 'vein_18_FirstMaskClean.png', isStep = true, cmap = cVein)

M0_Artery = setcmap(M0_ff_img, maskArtery, cmapArtery);
M0_Vein = setcmap(M0_ff_img, maskVein, cmapVein);
M0_AV = setcmap(M0_ff_img, maskArtery & maskVein, cmapAV);

M0_RGB = (M0_Artery + M0_Vein) .* ~(maskArtery & maskVein) + M0_AV + rescale(M0_ff_img) .* ~(maskArtery | maskVein);
saveImage(M0_RGB, ToolBox, 'all_19_RGB.png', isStep = true)

% 2)  Improvements of the first mask

if params.json.Mask.ImproveMask

    % 2) 0) Computation of the M0 in Diastole and in Systole

    [M0_Systole_img, M0_Diastole_img, M0_Systole_video] = compute_diasys(M0_video, maskArtery, 'mask');
    [M0_Sys_img, M0_Dia_img, ~] = compute_diasys(M0_ff_video, maskArtery, 'mask');
    saveImage(rescale(M0_Systole_img), ToolBox, 'artery_20_systole_img.png', isStep = true)
    saveImage(rescale(M0_Diastole_img), ToolBox, 'vein_20_diastole_img.png', isStep = true)

    % 2) 1) New Vesselness Mask

    Systole_Frangi = frangiVesselness(M0_Systole_img, 'artery_20', ToolBox);
    Diastole_Frangi = frangiVesselness(M0_Diastole_img, 'vein_20', ToolBox);
    Systole_Gabor = gaborVesselness(M0_Sys_img, ToolBox, 'artery_20');
    Diastole_Gabor = gaborVesselness(M0_Dia_img, ToolBox, 'vein_20');

    if params.json.Mask.VesselnessHolonet
        maskVesselness = maskVesselness & maskDiaphragm;
    else
        maskVesselness = (Systole_Frangi | Diastole_Frangi | Systole_Gabor | Diastole_Gabor) & maskDiaphragm;
    end

    maskVesselnessClean = removeDisconnected(maskVesselness, maskVesselness, maskCircle, 'all_20_VesselMask', ToolBox);

    % 2) 2) Diastole-Systole Image

    diasysArtery = M0_Systole_img - M0_Diastole_img;
    mDiasys = mean(diasysArtery, 'all', 'omitnan');
    diasysVein = mDiasys - diasysArtery;
    saveImage(diasysArtery, ToolBox, 'artery_21_diasys_img.png', isStep = true)
    saveImage(diasysVein, ToolBox, 'vein_21_diasys_img.png', isStep = true)

    RGBdiasys = labDuoImage(rescale(M0_Gabor), diasysArtery);
    saveImage(RGBdiasys, ToolBox, 'vessel_40_diasys_rgb.png', isStep = true)
    saveImage(RGBdiasys, ToolBox, 'DiaSysRGB.png')

    if diasysAnalysis % Systole/Diastole Analysis

        % 2) 3) Diastole-Systole based Segmentation
        [maskArtery, ~, quantizedImage] = processDiaSysSignal(diasysArtery, maskVesselnessClean, arteryParams, cArtery, 'artery_23');
        saveImage(maskArtery, ToolBox, 'artery_23_DiaSysMask.png', isStep = true, cmap = cArtery);
        saveImage(quantizedImage, ToolBox, 'artery_23_quantizedImage.png', isStep = true);
        quantizedImageRGB = quantizeImageToRGB(quantizedImage, arteryParams.classes);
        saveImage(quantizedImageRGB, ToolBox, 'artery_23_quantizedImage.png', isStep = true, cmap = cArtery);

        [~, maskVein, quantizedImage] = processDiaSysSignal(diasysVein, maskVesselnessClean, veinParams, cVein, 'vein_23');
        saveImage(maskVein, ToolBox, 'vein_23_DiaSysMask.png', isStep = true, cmap = cVein);
        saveImage(quantizedImage, ToolBox, 'vein_23_quantizedImage.png', isStep = true);
        quantizedImageRGB = quantizeImageToRGB(quantizedImage, veinParams.classes);
        saveImage(quantizedImageRGB, ToolBox, 'vein_23_quantizedImage.png', isStep = true, cmap = cVein);

    else % Second Correlation Analysis

        % 2) 3) Artery-Vein correlation based Segmentation
        [maskArtery, ~, R, ~, quantizedImage, level, color] = correlationSegmentation(M0_Systole_video, maskVesselnessClean, arteryParams);
        saveImage(R, ToolBox, 'artery_23_Correlation.png', isStep = true)
        RGBcorr = labDuoImage(M0_ff_img, R);
        saveImage(RGBcorr, ToolBox, 'artery_23_Correlation_rgb.png', isStep = true)
        graphThreshHistogram(R, level, maskVesselnessClean, color, 'artery_23');
        saveImage(quantizedImage, ToolBox, 'artery_23_quantizedImage.png', isStep = true);
        quantizedImageRGB = quantizeImageToRGB(quantizedImage, arteryParams.classes);
        saveImage(quantizedImageRGB, ToolBox, 'artery_23_quantizedImage_rgb.png', isStep = true, cmap = color);

        [~, maskVein, quantizedImage] = processDiaSysSignal(diasysVein, maskVesselnessClean, veinParams, cVein, 'vein_23');
        saveImage(maskVein, ToolBox, 'vein_23_DiaSysMask.png', isStep = true, cmap = cVein);
        saveImage(quantizedImage, ToolBox, 'vein_23_quantizedImage.png', isStep = true);
        quantizedImageRGB = quantizeImageToRGB(quantizedImage, veinParams.classes);
        saveImage(quantizedImageRGB, ToolBox, 'vein_23_quantizedImage_rgb.png', isStep = true);

    end

    % 3) Mask Clearing

    % 3) 0) Morphological Operations
    results = cell(2, 1);

    parfor i = 1:2

        if i == 1
            % Process artery mask
            results{i} = clearMasks(maskArtery, 'artery_30', cArtery, ToolBox);
        else
            % Process vein mask
            results{i} = clearMasks(maskVein, 'vein_30', cVein, ToolBox);
        end

    end

    maskArtery = results{1};
    maskVein = results{2};

    % 3) 1) Final Blob removal
    maskVessel = maskArtery | maskVein;
    maskArtery = removeDisconnected(maskArtery, maskVessel, maskCircle, 'artery_31_VesselMask', ToolBox);
    maskVein = removeDisconnected(maskVein, maskVessel, maskCircle, 'vein_31_VesselMask', ToolBox);

    % 3) 2) Force Create Masks in case they exist
    if (params.json.Mask.ForcedMasks == -1 || params.json.Mask.ForcedMasks == 1)

        if isfile(fullfile(ToolBox.path_main, 'mask', 'forceMaskArtery.png'))
            maskArtery = mat2gray(mean(imread(fullfile(ToolBox.path_main, 'mask', 'forceMaskArtery.png')), 3)) > 0;

            if size(maskArtery, 1) ~= maskCircle
                maskArtery = imresize(maskArtery, [numX, numY], "nearest");
            end

        elseif params.json.Mask.ForcedMasks == 1
            error("Cannot force Artery Mask because none given in the mask folder. Please create a forceMaskArtery.png file in the mask folder. (SET ForcedMasks to -1 to skip use the auto mask)");
        end

        if isfile(fullfile(ToolBox.path_main, 'mask', 'forceMaskVein.png'))
            maskVein = mat2gray(mean(imread(fullfile(ToolBox.path_main, 'mask', 'forceMaskVein.png')), 3)) > 0;

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

    % 3) 0) Morphological Operations
    results = cell(2, 1);

    parfor i = 1:2

        if i == 1
            % Process artery mask
            results{i} = clearMasks(maskArtery, 'artery_30', cArtery, ToolBox);
        else
            % Process vein mask
            results{i} = clearMasks(maskVein, 'vein_30', cVein, ToolBox);
        end

    end

    maskArtery = results{1};
    maskVein = results{2};

    maskVessel = maskArtery | maskVein;
    maskArtery = removeDisconnected(maskArtery, maskVessel, maskCircle, 'artery_31_VesselMask', ToolBox);
    maskVein = removeDisconnected(maskVein, maskVessel, maskCircle, 'vein_31_VesselMask', ToolBox);

end

% 3) 5) Create Vessel and Background Mask
maskArtery = maskArtery & maskDiaphragm;
maskVein = maskVein & maskDiaphragm;
maskVessel = maskArtery | maskVein;
maskBackground = not(maskVessel);

% 4) FINAL FIGURES

% 4) 1) RGB Figures
if isfile(fullfile(ToolBox.path_main, 'mask', 'forceMaskArtery.png')) || isfile(fullfile(ToolBox.path_main, 'mask', 'forceMaskVein.png'))

    maskArteryNoImport = results{1};
    maskVeinNoImport = results{2};

    maskVesselNoImport = maskArteryNoImport | maskVeinNoImport;
    maskArteryNoImport = removeDisconnected(maskArteryNoImport, maskVesselNoImport, maskCircle, 'artery_31_VesselMask', ToolBox);
    maskVeinNoImport = removeDisconnected(maskVeinNoImport, maskVesselNoImport, maskCircle, 'vein_31_VesselMask', ToolBox);

    M0_Artery = setcmap(M0_ff_img, maskArteryNoImport, cmapArtery);
    M0_Vein = setcmap(M0_ff_img, maskVeinNoImport, cmapVein);
    M0_AV = setcmap(M0_ff_img, maskArteryNoImport & maskVeinNoImport, cmapAV);
    M0_RGB = (M0_Artery + M0_Vein) .* ~(maskArteryNoImport & maskVeinNoImport) + M0_AV + rescale(M0_ff_img) .* ~(maskArteryNoImport | maskVeinNoImport);
    saveImage(M0_RGB, ToolBox, 'vessel_40_RGB_no_import.png', isStep = true)
end

M0_Artery = setcmap(M0_ff_img, maskArtery, cmapArtery);
M0_Vein = setcmap(M0_ff_img, maskVein, cmapVein);
M0_AV = setcmap(M0_ff_img, maskArtery & maskVein, cmapAV);

M0_RGB = (M0_Artery + M0_Vein) .* ~(maskArtery & maskVein) + M0_AV + rescale(M0_ff_img) .* ~(maskArtery | maskVein);
saveImage(M0_RGB, ToolBox, 'vessel_40_RGB.png', isStep = true)
saveImage(M0_RGB, ToolBox, 'RGB_img.png')

% 4) 2) Neighbours Mask

if params.json.Mask.AllNonVesselsAsBackground
    maskNeighbors = (maskBackground & ~maskVesselness) & maskDiaphragm;
else
    maskNeighbors = imdilate(maskArtery | maskVein, strel('disk', bgWidth)) & ~(maskArtery | maskVein);
end

cmapNeighbors = cmapLAB(256, [0 1 0], 0, [1 1 1], 1);

M0_Neighbors = setcmap(M0_ff_img, maskNeighbors, cmapNeighbors);

neighborsMaskSeg = (M0_Artery + M0_Vein) .* ~(maskArtery & maskVein) + ...
    M0_AV + M0_Neighbors + ...
    rescale(M0_ff_img) .* ~(maskArtery | maskVein | maskNeighbors);
saveImage(neighborsMaskSeg, ToolBox, 'neighbors_img.png')

% 4) 4) Save all images

saveImage(maskArtery, ToolBox, 'maskArtery.png')
saveImage(maskVein, ToolBox, 'maskVein.png')
saveImage(maskVessel, ToolBox, 'maskVessel.png')
saveImage(maskNeighbors, ToolBox, 'maskNeighbors.png')
saveImage(maskBackground, ToolBox, 'maskBackground.png')
saveImage(bwskel(maskArtery), ToolBox, 'skeletonArtery.png')
saveImage(bwskel(maskVein), ToolBox, 'skeletonVein.png')

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
