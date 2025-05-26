function [quantizedImageRGB] = quantizeImageToRGB(quantizedImage, classes)
% Function to convert a quantized image into an RGB image based on class labels
% Input:
%   quantizedImage: A 2D matrix where each pixel value corresponds to a class label
%   classes: A vector containing the unique class labels (e.g., [1, -1, 2, 0])
% Output:
%   quantizedImageRGB: A 3D RGB image where each class is represented by a specific color

numClasses = length(classes);
[numX, numY] = size(quantizedImage);

% Define RGB colors for each class
% Assuming classes are defined as follows:
% 1: Artery, -1: Vein, 2: Choroid, 0: Background
cArtery = [255 22 18] / 255;
cVein = [18 23 255] / 255;
cChoroid = [0 179 0] / 255;

% Create a mask for each class in the quantized image
% and assign the corresponding RGB color
quantizedRGB_Mask = repmat(quantizedImage, [1 1 3]);
quantizedImageRGB = zeros(numX, numY, 3);

% Assign colors based on the class labels
for i = 1:numClasses
    IDX = reshape([i + 1; i + 1; i + 1], 1, 1, 3);

    if classes(i) == 1
        quantizedImageRGB = quantizedImageRGB + reshape(cArtery, 1, 1, 3) .* (quantizedRGB_Mask == IDX);
    elseif classes(i) == -1
        quantizedImageRGB = quantizedImageRGB + reshape(cVein, 1, 1, 3) .* (quantizedRGB_Mask == IDX);
    elseif classes(i) == 2
        quantizedImageRGB = quantizedImageRGB + reshape(cChoroid, 1, 1, 3) .* (quantizedRGB_Mask == IDX);
    elseif classes(i) == 0
        quantizedImageRGB = quantizedImageRGB + reshape([1 1 1], 1, 1, 3) .* (quantizedRGB_Mask == IDX);
    end

end

end
