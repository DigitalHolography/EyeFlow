function createMasks(M0_ff, VesselSegmentationNet, AVSegmentationNet)
% createMasks - Creates masks for arteries, veins, and neighbors from a video of retinal images.
% Inputs:
%   M0_ff: 3D matrix of the video data (height x width x time)
% Outputs: (inside the ToolBox Cache)
%   maskArtery: Binary mask for arteries
%   maskVein: Binary mask for veins
%   maskNeighbors: Binary mask for neighboring regions

arguments
    M0_ff double
    VesselSegmentationNet = []
    AVSegmentationNet = []
end

ToolBox = getGlobalToolBox;
params = ToolBox.getParams;
saveFigures = params.saveFigures;
path = ToolBox.path_main;
folder_steps = fullfile('mask', 'steps');

if ~exist(fullfile(ToolBox.path_main, folder_steps), 'dir')
    mkdir(ToolBox.path_main, folder_steps);
end

% 0) Initialisation
[numX, numY, numFrames] = size(M0_ff);
xy_barycenter = ToolBox.Cache.xy_barycenter;
x_c = xy_barycenter(1);
y_c = xy_barycenter(2);

% Load mask parameters

if ~isfield(params.json, 'Mask')
    error('Mask parameters are missing in the JSON file.');
end

mask_params = params.json.Mask;
diaphragmRadius = mask_params.DiaphragmRadius;
cropChoroidRadius = mask_params.CropChoroidRadius;
vesselParams.threshold = mask_params.VascularThreshold;
vesselParams.classes = mask_params.VascularClasses;
forceVesselWidth = mask_params.ForceVesselWidth;
vesselnessMethod = mask_params.VesselSegmentationMethod;

% Load vesselness parameters
bgWidth = params.json.PulseAnalysis.LocalBackgroundWidth;
r1 = params.json.SizeOfField.SmallRadiusRatio;
r2 = params.json.SizeOfField.BigRadiusRatio;

cArtery = [0, 0, 0; 255, 22, 18] / 255;
cVein = [0, 0, 0; 18, 23, 255] / 255;

cmapArtery = ToolBox.Cache.cmapArtery;
cmapVein = ToolBox.Cache.cmapVein;
cmapAV = ToolBox.Cache.cmapAV;

M0_ff_img = squeeze(mean(M0_ff, 3));
saveMaskImage(M0_ff_img, 'all_10_M0.png', isStep = true)

maskDiaphragm = diskMask(numX, numY, diaphragmRadius);
saveMaskImage(rescale(M0_ff_img) + maskDiaphragm .* 0.5, 'all_11_maskDiaphragm.png', isStep = true)
maskCircle = diskMask(numX, numY, cropChoroidRadius, 'center', [x_c / numX, y_c / numY]);

scoreMaskArtery = NaN;
scoreMaskVein = NaN;

if mask_params.AutoCompute

    % 1) First Masks and Correlation
    M0_video = M0_ff;
    A = ones(1, 1, numFrames);
    B = A .* maskDiaphragm;
    M0_video(~B) = NaN;
    clear A B

    M0_img = squeeze(mean(M0_video, 3, 'omitnan'));

    % 1) 1) Compute vesselness response

    switch vesselnessMethod

        case 'matchedFilter'
            fprintf("Compute vesselness using matched filter\n");

            % Matched Vesselness
            [~, maskVesselness] = matchedFilterVesselDetection(M0_img, ...
                'threshold', 0.6);

            saveMaskImage(maskVesselness, 'all_12_matched_filter_mask.png', isStep = true)

        case 'frangi'
            fprintf("Compute vesselness using Frangi filter\n");
            % Frangi Vesselness
            [maskVesselness, M0_Frangi] = frangiVesselness(M0_img, ...
                'range', [4, 6], 'step', 1);

            saveMaskImage(maskVesselness, 'all_12_frangi_mask.png', isStep = true)
            saveMaskImage(M0_Frangi, 'all_12_frangi_img.png', isStep = true)

        case 'gabor'
            fprintf("Compute vesselness using Gabor filter\n");
            % Gabor Vesselness
            [maskVesselness, M0_Gabor] = gaborVesselness(M0_ff_img, ...
                'range', [4, 6], 'step', 1);

            saveMaskImage(maskVesselness & maskDiaphragm, 'all_12_gabor_mask.png', isStep = true)
            saveMaskImage(M0_Gabor, 'all_12_gabor_img.png', isStep = true)

        case 'AI'
            fprintf("Compute vesselness using SegmentationNet\n");
            % SegmentationNet Vesselness
            maskVesselness = getSegmentationNetVesselness(M0_ff_img, VesselSegmentationNet);
            saveMaskImage(maskVesselness, 'all_14_maskSegmentationNet.png', isStep = true)

    end

    % 1) 2) Clean vesselness response
    maskVesselnessClean = maskVesselness & bwareafilt(maskVesselness | maskCircle, 1, 8) & maskDiaphragm;
    saveMaskImage(maskVesselnessClean + maskCircle * 0.5, 'all_15_VesselMask_clear.png', isStep = true)

    % 1) 3) Pre-mask arteries using intensity information
    [maskArteryTmp, maskVeinTmp] = preMaskArtery(M0_ff, maskVesselnessClean);
    saveMaskImage(maskArteryTmp, 'artery_16_PreMask.png', isStep = true, cmap = cArtery);
    saveMaskImage(maskVeinTmp, 'vein_16_PreMask.png', isStep = true, cmap = cVein);

    % 1) 4) If using SegmentationNet, compute correlation and/or dia/sys to obtain artery vein masks

    if mask_params.AVCorrelationSegmentationNet || mask_params.AVDiasysSegmentationNet

        % Compute artery/vein masks using SegmentationNet
        [maskArtery, maskVein, scoreMaskArtery, scoreMaskVein] = createMasksSegmentationNet(M0_ff, M0_ff_img, AVSegmentationNet, maskArteryTmp);
        saveMaskImage(maskVein, 'vein_21_SegmentationNet.png', isStep = true, cmap = cVein);
        saveMaskImage(maskArtery, 'artery_21_SegmentationNet.png', isStep = true, cmap = cVein);

    else
        % Compute artery/vein masks using correlation-based segmentation
        [maskArtery, maskVein, R_ArterialSignal, arterialSignal, quantizedImage, level, color] = ...
            correlationSegmentation(M0_video, maskArteryTmp, maskVesselnessClean, vesselParams);

        % 1) 5) Save all images
        if saveFigures
            t = ToolBox.Cache.t;
            graphSignal('all_15_arterialSignal', t, squeeze(arterialSignal), '-', color(end, :), ...
                Title = 'Arterial Signal', xlabel = 'Time(s)', ylabel = 'Power Doppler (a.u.)');

            saveMaskImage(R_ArterialSignal, 'all_15_Correlation.png', isStep = true)
            RGBcorr = labDuoImage(M0_ff_img, R_ArterialSignal);
            saveMaskImage(RGBcorr, 'all_15_Correlation_rgb.png', isStep = true)

            [quantizedImageRGB] = quantizeImageToRGB(quantizedImage, vesselParams.classes);
            saveMaskImage(quantizedImageRGB, 'all_16_quantizedImage.png', isStep = true, cmap = color);
            graphThreshHistogram(R_ArterialSignal, level, maskVesselnessClean, color, 'all_16')
        end

    end

    % 3) Mask Clearing

    % 3) 0) Morphological Operations

    % Ensure masks are within diaphragm
    maskArtery = maskArtery & maskDiaphragm;
    maskVein = maskVein & maskDiaphragm;

    [mask_dilated, mask_closed, mask_opened, mask_widened] = clearMasks(maskArtery, ...
        'min_area', mask_params.MinPixelSize, ...
        'imclose_radius', mask_params.ImcloseRadius, ...
        'min_width', mask_params.MinimumVesselWidth, ...
        'imdilate_size', mask_params.FinalDilation);

    saveMaskImage(mask_dilated, 'artery_20_ClearedMask.png', isStep = true, cmap = cArtery);
    saveMaskImage(mask_opened, 'artery_20_ClearedMask_opened.png', isStep = true, cmap = cArtery);
    saveMaskImage(mask_closed, 'artery_20_ClearedMask_closed.png', isStep = true, cmap = cArtery);
    saveMaskImage(mask_widened, 'artery_20_ClearedMask_widened.png', isStep = true, cmap = cArtery);

    maskArtery = mask_dilated;

    [mask_dilated, mask_closed, mask_opened, mask_widened] = clearMasks(maskVein, ...
        'min_area', mask_params.MinPixelSize, ...
        'imclose_radius', mask_params.ImcloseRadius, ...
        'min_width', mask_params.MinimumVesselWidth, ...
        'imdilate_size', mask_params.FinalDilation);

    saveMaskImage(mask_dilated, 'vein_20_ClearedMask.png', isStep = true, cmap = cVein);
    saveMaskImage(mask_opened, 'vein_20_ClearedMask_opened.png', isStep = true, cmap = cVein);
    saveMaskImage(mask_closed, 'vein_20_ClearedMask_closed.png', isStep = true, cmap = cVein);
    saveMaskImage(mask_widened, 'vein_20_ClearedMask_widened.png', isStep = true, cmap = cVein);

    maskVein = mask_dilated;

    % 3) 1) Final Blob removal
    maskVessel = maskArtery | maskVein;

    maskArtery = maskArtery & bwareafilt(maskVessel | maskCircle, 1, 8);
    saveMaskImage(maskArtery + maskCircle * 0.5, 'artery_21_VesselMask_clear.png', isStep = true)

    maskVein = maskVein & bwareafilt(maskVessel | maskCircle, 1, 8);
    saveMaskImage(maskVein + maskCircle * 0.5, 'vein_21_VesselMask_clear.png', isStep = true)

end

% 3) 2) a) Look for a target mask to register from
if (mask_params.RegisteredMasks == -1 || mask_params.RegisteredMasks == 1)

    if (isfile(fullfile(path, 'mask', 'similarMaskArtery.png')) ...
            && isfile(fullfile(path, 'mask', 'similarM0.png'))) ...
            || (isfile(fullfile(path, 'mask', 'similarMaskVein.png')) ...
            && isfile(fullfile(path, 'mask', 'similarM0.png')))
        M0_ff_img = squeeze(mean(M0_ff, 3));
        similarM0 = mat2gray(mean(imread(fullfile(path, 'mask', 'similarM0.png')), 3));
        [ux, uy] = nonrigidregistration(similarM0, M0_ff_img, fullfile(ToolBox.path_png, folder_steps), 'Reg');

    end

    if isfile(fullfile(path, 'mask', 'similarMaskArtery.png')) ...
            && isfile(fullfile(path, 'mask', 'similarM0.png'))
        similarMaskArtery = mat2gray(mean(imread(fullfile(path, 'mask', 'similarMaskArtery.png')), 3)) > 0;
        similarMaskArtery = imresize(similarMaskArtery, [numX, numY], "nearest");
        [Xq, Yq] = meshgrid(1:numX, 1:numY);
        maskArtery = interp2(single(similarMaskArtery), Xq + ux, Yq + uy, 'linear', 0) > 0;
    end

    if isfile(fullfile(path, 'mask', 'similarMaskVein.png')) ...
            && isfile(fullfile(path, 'mask', 'similarM0.png'))
        similarMaskVein = mat2gray(mean(imread(fullfile(path, 'mask', 'similarMaskVein.png')), 3)) > 0;
        similarMaskVein = imresize(similarMaskVein, [numX, numY], "nearest");
        [Xq, Yq] = meshgrid(1:numX, 1:numY);
        maskVein = interp2(single(similarMaskVein), Xq + ux, Yq + uy, 'linear', 0) > 0;
    end

end

% 3) 2) b) Force Create Masks in case they exist
if (mask_params.ForcedMasks == -1 || mask_params.ForcedMasks == 1)

    if isfile(fullfile(path, 'mask', 'forceMaskArtery.png'))
        maskArtery = mat2gray(mean(imread(fullfile(path, 'mask', 'forceMaskArtery.png')), 3)) > 0;

        if size(maskArtery, 1) ~= maskCircle
            maskArtery = imresize(maskArtery, [numX, numY], "nearest");
        end

    elseif mask_params.ForcedMasks == 1
        error("Cannot force Artery Mask because none given in the mask folder. Please create a forceMaskArtery.png file in the mask folder. (SET ForcedMasks to -1 to skip use the auto mask)");
    end

    if isfile(fullfile(path, 'mask', 'forceMaskVein.png'))
        maskVein = mat2gray(mean(imread(fullfile(path, 'mask', 'forceMaskVein.png')), 3)) > 0;

        if size(maskVein, 1) ~= maskCircle
            maskVein = imresize(maskVein, [numX, numY], "nearest");
        end

    elseif mask_params.ForcedMasks == 1
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

% 3) 5) Create Vessel and Background Mask
maskArtery = maskArtery & maskDiaphragm;
maskVein = maskVein & maskDiaphragm;
maskVessel = maskArtery | maskVein;
maskAV = maskArtery & maskVein;
maskBackground = not(maskVessel);

% 3) 6) Neighbours Mask

bgDist = params.json.PulseAnalysis.LocalBackgroundDist;

maskVessel_tmp = imdilate(maskArtery | maskVein, strel('disk', bgDist));
maskNeighbors = imdilate(maskVessel_tmp, strel('disk', bgWidth)) & ~(maskVessel_tmp);

maskNeighbors = maskNeighbors & maskDiaphragm;

% 4) FINAL FIGURES

% Normalize and prepare base image
M0_ff_img = rescale(M0_ff_img);

% Apply color maps to each region
M0_Artery = setcmap(M0_ff_img, maskArtery, cmapArtery);
M0_Vein = setcmap(M0_ff_img, maskVein, cmapVein);
M0_AV = setcmap(M0_ff_img, maskAV, cmapAV);
M0_bkg = M0_ff_img .* maskBackground;

% Combine into final RGB image
M0_RGB = M0_Artery + M0_Vein + M0_AV + M0_bkg;

% Save main mask images
saveMaskImage(M0_Artery + M0_ff_img .* ~maskArtery, 'M0_Artery.png', ForceFigure = true)
saveMaskImage(M0_Vein + M0_ff_img .* ~maskVein, 'M0_Vein.png', ForceFigure = true)
saveMaskImage(M0_RGB, 'M0_RGB.png', ForceFigure = true)

if saveFigures

    % Neighbors Mask Figures
    cmapNeighbors = cmapLAB(256, [0 1 0], 0, [1 1 1], 1);
    M0_Neighbors = setcmap(M0_ff_img, maskNeighbors, cmapNeighbors);

    neighborsMaskSeg = M0_Artery + M0_Vein + M0_AV + M0_Neighbors + ...
        M0_ff_img .* ~(maskArtery | maskVein | maskNeighbors);

    saveMaskImage(neighborsMaskSeg, 'neighbors_img.png')

    % Save all images

    saveMaskImage(maskArtery, 'maskArtery.png')
    saveMaskImage(maskVein, 'maskVein.png')
    saveMaskImage(maskVessel, 'maskVessel.png')
    saveMaskImage(maskNeighbors, 'maskNeighbors.png')
    saveMaskImage(maskBackground, 'maskBackground.png')
    saveMaskImage(bwskel(maskArtery), 'skeletonArtery.png')
    saveMaskImage(bwskel(maskVein), 'skeletonVein.png')

    % Mask Section & Force Barycenter
    createMaskSection(ToolBox, M0_ff_img, r1, r2, xy_barycenter, 'vessel_map_artery', maskArtery, thin = 0.01);
    createMaskSection(ToolBox, M0_ff_img, r1, r2, xy_barycenter, 'vessel_map_vein', [], maskVein, thin = 0.01);
    createMaskSection(ToolBox, M0_ff_img, r1, r2, xy_barycenter, 'vessel_map', maskArtery, maskVein, thin = 0.01);

end

% 5) Save masks in Cache
ToolBox.Cache.maskArtery = maskArtery;
ToolBox.Cache.maskVein = maskVein;
ToolBox.Cache.maskNeighbors = maskNeighbors;
ToolBox.Cache.maskBackground = maskBackground;
ToolBox.Cache.xy_barycenter = xy_barycenter;
ToolBox.Cache.M0_RGB = M0_RGB;
ToolBox.Cache.scoreMaskArtery = scoreMaskArtery;
ToolBox.Cache.scoreMaskVein = scoreMaskVein;

close all
end
