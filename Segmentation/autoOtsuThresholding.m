function [quantizedImage, level, color] = autoOtsuThresholding(image, mask, classes)
% autoOtsuThresholding - Segments an image using Otsu's method and returns quantized image, thresholds, and colors.
% Inputs:
% image = a 2D matrix where the segmentation is performed
% mask = a binary mask where the segmentation is performed
% classes = [-1 1] a 0 and 1 matrix where 1 is the class selected and 0 the rejected
% Outputs:
% quantizedImage = a 2D matrix where the segmentation is performed
% level = the thresholds used for the segmentation
% color = the colors used for the segmentation

if nargin < 3
    classes = [-1 1]; % Default classes: 0 = background, 1 = artery, -1 = vein, 2 = choroid
end

if nargin < 2
    mask = true(size(image)); % Default mask is the whole image
end

cArtery = [255 22 18] / 255;
cVein = [18 23 255] / 255;
cChoroid = [0 179 0] / 255;

numClasses = length(classes);
color = zeros(numClasses, 3);

for i = 1:numClasses

    switch classes(i)
        case 1
            color(i, :) = cArtery;
        case -1
            color(i, :) = cVein;
        case 0
            color(i, :) = [1 1 1];
        case 2
            color(i, :) = cChoroid;
    end

end

image = rescale(image);
image(isnan(image)) = 0;
level = multithresh(image(mask), numClasses - 1);
quantizedImage = imquantize(image - 2 * ~mask, [-1 level]);

end
