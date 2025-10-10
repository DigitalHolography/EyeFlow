function [maskArtery, maskVein, R, signalArtery, quantizedImage, level, color] = correlationSegmentation(video, preMaskArtery, maskVesselness, params)
% correlationSegmentation - Segments vessels from a video using the provided parameters.
% Inputs:
%   video: 3D matrix of the video data (height x width x time)
%   mask: 2D binary mask indicating vesselness (height x width)
%   params: structure containing parameters for segmentation
%       .threshold: manual threshold for segmentation (-1 to 1)
%       .classes: number of classes for automatic Otsu segmentation
% Outputs:
%   maskArtery: binary mask for arteries
%   maskVein: binary mask for veins
%   R: correlation map of arterial signal
%   signalArtery: 1D signal representing the arterial activity
%   quantizedImage: quantized image after Otsu segmentation

% Initialize output variables
classes = params.classes;
numClasses = length(classes);

%  1) Compute first correlation
% compute signal in 3 dimentions for correlation in the mask
signalArtery = sum(video .* preMaskArtery, [1 2], 'omitnan');
signalArtery = signalArtery ./ nnz(preMaskArtery);

outlier_frames_mask = isoutlier(signalArtery, "movmedian", 5, "ThresholdFactor", 2);
video = interpolateOutlierFrames(video, outlier_frames_mask);

signalArtery = sum(video .* preMaskArtery, [1 2], 'omitnan');
signalArtery = signalArtery ./ nnz(preMaskArtery);

% compute local-to-average signal wave zero-lag correlation
signal_centered = signalArtery - mean(signalArtery, 3, 'omitnan');
video_centered = video - mean(video, 'all', 'omitnan');
R = mean(video_centered .* signal_centered, 3) ./ (std((video_centered), [], 'all', 'omitnan') * std(signal_centered, [], 3));

% 2) Segment Vessels
if isempty(params.threshold)
    % Automatic Otsu segmentation is performed
    [quantizedImage, level, color] = autoOtsuThresholding(R, maskVesselness, params.classes);

    % Initialize masks
    maskArtery = false(size(R));
    maskVein = false(size(R));

    % Process each class
    for i = 1:numClasses

        if classes(i) == 1
            maskArtery = maskArtery | (quantizedImage == i + 1);
        elseif classes(i) == -1
            maskVein = maskVein | (quantizedImage == i + 1);
        end

    end

elseif abs(params.threshold) <= 1

    % Create masks based on threshold
    maskArtery = (R > params.threshold) & maskVesselness;
    maskVein = (R < params.threshold) & maskVesselness;

    % Create quantized image
    quantizedImage = zeros(size(R));
    quantizedImage(maskVein) = 0.5; % Veins
    quantizedImage(maskArtery) = 1; % Arteries

    % Assign colors for visualization
    color = [255, 22, 18; % Red for arteries
             18, 23, 255] ./ 255; % Blue for veins
    level = [params.threshold]; % Threshold levels for visualization
else

    % Automatic Otsu segmentation is performed
    [quantizedImage, level, color] = autoOtsuThresholding(R, maskVesselness, params.classes);

    % Initialize masks
    maskArtery = false(size(R));
    maskVein = false(size(R));

    % Process each class
    for i = 1:numClasses

        if classes(i) == 1
            maskArtery = maskArtery | (quantizedImage == i + 1);
        elseif classes(i) == -1
            maskVein = maskVein | (quantizedImage == i + 1);
        end

    end

end

end
