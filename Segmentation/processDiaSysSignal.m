function [mask1, mask2, quantizedImage] = processDiaSysSignal(diasys, maskClean, params, cmap, suffix)
% This function processes the Diastolic-Systolic signal to segment arteries
% and veins based on a provided threshold or using Otsu's method.
% Inputs:
%   diasys: 2D matrix of the Diastolic-Systolic signal
%   maskClean: 2D binary mask indicating the area of interest
%   params: structure containing parameters for segmentation
%       .threshold: manual threshold for segmentation (-1 to 1)
%       .classes: number of classes for automatic Otsu segmentation
%   cmap: colormap for visualization
%   suffix: string to append to the output files
% Outputs:
%   mask1: binary mask for arteries
%   mask2: binary mask for veins

classes = params.classes;
numClasses = length(classes);

if isempty(params.threshold)

    % Automatic Otsu thresholding
    [quantizedImage, level, color] = autoOtsuThresholding(diasys, maskClean, params.classes);

    % Initialize masks
    mask1 = false(size(diasys));
    mask2 = false(size(diasys));

    % Process each class
    for i = 1:numClasses

        if classes(i) == 1
            mask1 = mask1 | (quantizedImage == i + 1); % Arteries
        elseif classes(i) == -1
            mask2 = mask2 | (quantizedImage == i + 1); % Veins
        end

    end

    % Visualize results
    graphThreshHistogram(diasys, level, maskClean, color, suffix);

elseif abs(params.threshold) <= 1

    % Rescale input data to [0, 1] range
    diasys = rescale(diasys);

    % Create masks based on threshold
    mask1 = diasys >= params.threshold;
    mask2 = diasys <= params.threshold;

    % Create quantized image
    quantizedImage = zeros(size(diasys));
    quantizedImage(mask1) = 1; % Arteries
    quantizedImage(mask2) = 1/2; % Veins

    % Visualize results
    graphThreshHistogram(diasys, params.threshold, maskClean, cmap, suffix);

else

    % Automatic Otsu thresholding
    [quantizedImage, level, color] = autoOtsuThresholding(diasys, maskClean, params.classes);

    % Initialize masks
    mask1 = false(size(diasys));
    mask2 = false(size(diasys));

    % Process each class
    for i = 1:numClasses

        if classes(i) == 1
            mask1 = mask1 | (quantizedImage == i + 1); % Arteries
        elseif classes(i) == -1
            mask2 = mask2 | (quantizedImage == i + 1); % Veins
        end

    end

    % Visualize results
    graphThreshHistogram(diasys, level, maskClean, color, suffix);

end

end
