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

if ~isempty(params.threshold)
    % Manual threshold
    diasys = rescale(diasys);
    mask1 = diasys >= params.threshold;
    mask2 = diasys <= params.threshold;
    graphThreshHistogram(diasys, params.threshold, maskClean, cmap, suffix);
    quantizedImage = zeros(size(diasys));
    quantizedImage(mask1) = 1/2; % Arteries
    quantizedImage(mask2) = 1; % Veins
else
    % Automatic Otsu thresholding
    [quantizedImage, level, color] = autoOtsuThresholding(diasys, maskClean, params.classes);

    classes = params.classes;
    numClasses = length(classes);
    [numX, numY] = size(diasys);

    mask1 = zeros(numX, numY, 'logical');
    mask2 = zeros(numX, numY, 'logical');

    for i = 1:numClasses

        if classes(i) == 1
            mask1 = mask1 | (quantizedImage == i + 1);
        elseif classes(i) == -1
            mask2 = mask2 | (quantizedImage == i + 1);
        end

    end

    graphThreshHistogram(diasys, level, maskClean, color, suffix);

end

end
