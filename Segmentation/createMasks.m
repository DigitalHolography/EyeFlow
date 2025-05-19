function [maskArtery, maskVein, maskNeighbors] = createMasks(M0_ff_video, xy_barycenter)

ToolBox = getGlobalToolBox;
params = ToolBox.getParams;

if ~isfolder(fullfile(ToolBox.path_png, 'mask'))
    mkdir(ToolBox.path_png, 'mask')
    mkdir(ToolBox.path_eps, 'mask')
    mkdir(fullfile(ToolBox.path_png, 'mask'), 'steps')
    mkdir(fullfile(ToolBox.path_eps, 'mask'), 'steps')
end

folder_steps = fullfile('mask', 'steps');

%% 0) Initialisation

% 0) 1) Parameters Initialisation

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

%% 1) First Masks and Correlation

maskDiaphragm = diskMask(numX, numY, diaphragmRadius);

M0_ff_img = squeeze(mean(M0_ff_video, 3));
M0_ff_video_centered = M0_ff_video - mean(M0_ff_video, [1 2]);
saveImage(M0_ff_img, ToolBox, 'all_10_M0.png', isStep = true)

if ~isfile(fullfile(ToolBox.path_gif, sprintf("%s_M0.gif", ToolBox.folder_name)))
    writeGifOnDisc(imresize(rescale(M0_ff_video), 0.5), "M0")
end

saveImage(rescale(M0_ff_img) + maskDiaphragm .* 0.5, ToolBox, 'all_11_maskDiaphragm.png', isStep = true)

% 1) 1) Compute vesselness response

[maskVesselnessFrangi] = frangiVesselness(M0_ff_img, 'all_12', ToolBox);
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
% compute signal in 3 dimentions for correlation in all vessels
vascularSignal = sum(M0_ff_video .* maskVesselnessClean, [1 2]);
vascularSignal = vascularSignal ./ nnz(maskVesselnessClean);

t = linspace(0, numFrames * ToolBox.stride / ToolBox.fs / 1000, numFrames);
tLabel = 'Time(s)';
yLabel = 'Power Doppler (a.u.)';

graphSignal('all_15_vascularSignal', folder_steps, t, squeeze(vascularSignal), '-', cVascular, Title = 'Vascular Signal', xlabel = tLabel, ylabel = yLabel);

% compute local-to-average signal wave zero-lag correlation
vascularSignal_centered = vascularSignal - mean(vascularSignal, 3);
R_VascularSignal = mean(M0_ff_video_centered .* vascularSignal_centered, 3) ./ (std((M0_ff_video_centered), [], 3) * std(vascularSignal_centered, [], 3));
saveImage(R_VascularSignal, ToolBox, 'all_15_Correlation.png', isStep = true)

mR_vascular = sum(R_VascularSignal .* fftshift(diskMask(numX, numY, 0.1)), [1 2]) ./ nnz(fftshift(diskMask(numX, numY, 0.1)));

RGBcorr = labDuoImage(M0_ff_img, R_VascularSignal - mR_vascular);
saveImage(RGBcorr, ToolBox, 'all_15_Correlation_rgb.png', isStep = true)

% 1) 4) Segment Vessels

cArtery = [255 22 18] / 255;
cVein = [18 23 255] / 255;

cmapArtery = [0 0 0; cArtery];
cmapVein = [0 0 0; cVein];
cmapVessels = [cVein; cArtery];

if vesselParams.threshold >= -1 && vesselParams.threshold <= 1
    % IF Manual Thresholds have been set between -1 and 1 then they are used

    maskArtery = (R_VascularSignal > vesselParams.threshold) .* maskVesselnessClean;
    maskVein = (R_VascularSignal < vesselParams.threshold) .* maskVesselnessClean;
    graphThreshHistogram(R_VascularSignal, vesselParams.threshold, maskVesselnessClean, cmapVessels, 'all_16')

else
    % ELSE automatic Otsu segmentation is performed
    % Number of classes for Vessels: 4
    % 1 & 2 = Veins & CoroidalVessels, 3 = CoroidalVessel, 4 = Arteries
    [maskArtery, maskVein] = autoOtsuThresholding(R_VascularSignal, maskVesselnessClean, vesselParams.classes, 'all_16');
end

saveImage(maskArtery, ToolBox, 'artery_17_FirstMask.png', isStep = true, cmap = cmapArtery)
saveImage(maskVein, ToolBox, 'vein_17_FirstMask.png', isStep = true, cmap = cmapVein)

% Remove small blobs
maskArtery = bwareaopen(maskArtery, minPixelSize);
maskVein = bwareaopen(maskVein, minPixelSize);
maskChoroid = maskVesselness & ~(maskArtery | maskVein);

saveImage(maskArtery, ToolBox, 'artery_18_FirstMaskClean.png', isStep = true, cmap = cmapArtery)
saveImage(maskVein, ToolBox, 'vein_18_FirstMaskClean.png', isStep = true, cmap = cmapVein)

RGBM0(:, :, 1) = rescale(M0_ff_img) + maskArtery;
RGBM0(:, :, 2) = rescale(M0_ff_img) + maskChoroid;
RGBM0(:, :, 3) = rescale(M0_ff_img) + maskVein;
saveImage(RGBM0, ToolBox, 'all_19_RGB.png', isStep = true)

%% 2)  Improvements of the first mask

if params.json.Mask.ImproveMask

    % 2) 0) Computation of the M0 in Diastole and in Systole

    [M0_Systole_img, M0_Diastole_img, M0_Systole_video] = compute_diasys(M0_ff_video, maskArtery, 'mask');
    saveImage(rescale(M0_Systole_img), ToolBox, 'artery_20_systole_img.png', isStep = true)
    saveImage(rescale(M0_Diastole_img), ToolBox, 'vein_20_diastole_img.png', isStep = true)

    % 2) 1) New Vesselness Mask

    Systole_Frangi = frangiVesselness(M0_Systole_img, 'artery_20', ToolBox);
    Diastole_Frangi = frangiVesselness(M0_Diastole_img, 'vein_20', ToolBox);
    Systole_Gabor = gaborVesselness(M0_Systole_img, ToolBox, 'artery_20');
    Diastole_Gabor = gaborVesselness(M0_Diastole_img, ToolBox, 'vein_20');

    if params.json.Mask.VesselnessHolonet
        maskVesselness = maskVesselness & maskDiaphragm;
    else
        maskVesselness = (Systole_Frangi | Diastole_Frangi | Systole_Gabor | Diastole_Gabor) & maskDiaphragm;
    end

    maskVesselnessClean = removeDisconnected(maskVesselness, maskVesselness, maskCircle, 'all_20_VesselMask', ToolBox);

    % 2) 2) Diastole-Systole Image

    diasysArtery = M0_Systole_img - M0_Diastole_img;
    mDiasys = sum(diasysArtery .* fftshift(diskMask(numX, numY, 0.1)), [1 2]) ./ nnz(fftshift(diskMask(numX, numY, 0.1)));
    diasysVein = mDiasys - diasysArtery;
    saveImage(diasysArtery, ToolBox, 'artery_21_diasys_img.png', isStep = true)
    saveImage(diasysVein, ToolBox, 'vein_21_diasys_img.png', isStep = true)

    RGBdiasys = labDuoImage(rescale(M0_Gabor), (diasysArtery - mDiasys));
    saveImage(RGBdiasys, ToolBox, 'vessel_40_diasys_rgb.png', isStep = true)
    saveImage(RGBdiasys, ToolBox, 'DiaSysRGB.png')

    if diasysAnalysis % Systole/Diastole Analysis

        % 2) 3) Diastole-Systole based Segmentation
        maskArtery = processDiaSysSignal(diasysArtery, maskVesselnessClean, arteryParams, cmapArtery, 'artery_23');
        [~, maskVein] = processDiaSysSignal(diasysVein, maskVesselnessClean, veinParams, cmapVein, 'vein_23');

    else % Second Correlation Analysis

        % 2) 3) Artery-Vein correlation based Segmentation
        maskArtery = processVascularSignal(M0_Systole_video, maskArtery, maskVesselnessClean, arteryParams, cmapArtery, 'artery_23', ToolBox);
        %         maskVein = processVascularSignal(M0_Diastole_video, maskVein, maskVesselnessClean, veinParams, cmapVein, 'vein_23', ToolBox);
        [~, maskVein] = processDiaSysSignal(diasysVein, maskVesselnessClean, veinParams, cmapVein, 'vein_23');

    end

    %% 3) Mask Clearing

    % 3) 0) Morphological Operations
    results = cell(2, 1);

    parfor i = 1:2

        if i == 1
            % Process artery mask
            results{i} = clearMasks(maskArtery, 'artery_30', cmapArtery, ToolBox);
        else
            % Process vein mask
            results{i} = clearMasks(maskVein, 'vein_30', cmapVein, ToolBox);
        end

    end

    maskArtery = results{1};
    maskVein = results{2};

    % 3) 1) Final Blob removal
    maskVessel = maskArtery | maskVein;
    maskArtery = removeDisconnected(maskArtery, maskVessel, maskCircle, 'artery_31_VesselMask', ToolBox);
    maskVein = removeDisconnected(maskVein, maskVessel, maskCircle, 'vein_31_VesselMask', ToolBox);

    % 3) 1 prime) HoloNet intervention
    if params.json.Mask.ChoroidHolonet

        try
            holonet_vessels = getHolonetprediction(M0_ff_img);
        catch
            warning("The Holonet ONNX-model couldn't be found.")
            holonet_vessels = maskVesselnessClean;
        end

        maskArtery = maskArtery & holonet_vessels;
        maskVein = maskVein & holonet_vessels;
    end

    % 3) 2) Force Create Masks in case they exist

    if isfile(fullfile(ToolBox.path_main, 'mask', 'forceMaskArtery.png'))
        maskArtery = mat2gray(mean(imread(fullfile(ToolBox.path_main, 'mask', 'forceMaskArtery.png')), 3)) > 0;

        if size(maskArtery, 1) ~= maskCircle
            maskArtery = imresize(maskArtery, [numX, numY], "nearest");
        end

    end

    if isfile(fullfile(ToolBox.path_main, 'mask', 'forceMaskVein.png'))
        maskVein = mat2gray(mean(imread(fullfile(ToolBox.path_main, 'mask', 'forceMaskVein.png')), 3)) > 0;

        if size(maskVein, 1) ~= maskCircle
            maskVein = imresize(maskVein, [numX, numY], "nearest");
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
            results{i} = clearMasks(maskArtery, 'artery_30', cmapArtery, ToolBox);
        else
            % Process vein mask
            results{i} = clearMasks(maskVein, 'vein_30', cmapVein, ToolBox);
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

%% 4) FINAL FIGURES

% 4) 1) RGB Figures
cmapArtery = cmapLAB(256, [0 0 0], 0, [1 0 0], 1/3, [1 1 0], 2/3, [1 1 1], 1);
cmapVein = cmapLAB(256, [0 0 0], 0, [0 0 1], 1/3, [0 1 1], 2/3, [1 1 1], 1);
cmapAV = cmapLAB(256, [0 0 0], 0, [1 0 1], 1/3, [1 1 1], 1);

M0_Artery = setcmap(M0_ff_img, maskArtery, cmapArtery);
M0_Vein = setcmap(M0_ff_img, maskVein, cmapVein);
M0_AV = setcmap(M0_ff_img, maskArtery & maskVein, cmapAV);

M0_RGB = (M0_Artery + M0_Vein) .* ~(maskArtery & maskVein) + M0_AV + rescale(M0_ff_img) .* ~(maskArtery | maskVein);
saveImage(M0_RGB, ToolBox, 'vessel_40_RGB.png', isStep = true)
saveImage(M0_RGB, ToolBox, 'RGB_img.png')

if isfile(fullfile(ToolBox.path_main, 'mask', 'forceMaskArtery.png')) || isfile(fullfile(ToolBox.path_main, 'mask', 'forceMaskVein.png'))
    M0_Artery = setcmap(M0_ff_img, results{1}, cmapArtery);
    M0_Vein = setcmap(M0_ff_img, results{2}, cmapVein);
    M0_AV = setcmap(M0_ff_img, results{1} & results{2}, cmapAV);
    M0_RGB = (M0_Artery + M0_Vein) .* ~(results{1} & results{2}) + M0_AV + rescale(M0_ff_img) .* ~(results{1} | results{2});
    saveImage(M0_RGB, ToolBox, 'vessel_40_RGB_no_import.png', isStep = true)
end

% 4) 2) Neighbours Mask

if params.json.Mask.AllNonVesselsAsBackground
    maskNeighbors = (maskBackground & ~maskVesselness) & maskDiaphragm;
else
    maskNeighbors = imdilate(maskArtery | maskVein, strel('disk', bgWidth)) - (maskArtery | maskVein);
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
try
    getLongestArteryBranch(maskArtery, xy_barycenter, 'Artery');
    getLongestArteryBranch(maskVein, xy_barycenter, 'Vein');
catch ME

    for i = 1:length(ME.stack)
        disp("Error in getLongestArteryBranch: " + ME.stack(i).name + " at line " + ME.stack(i).line)
    end

end

close all
end
