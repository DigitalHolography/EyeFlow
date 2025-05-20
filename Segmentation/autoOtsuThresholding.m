function [maskArtery, maskVein, maskChoroid, maskBackground] = autoOtsuThresholding(image, mask, classes, name)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% classes = [0 0 1 1] a 0 and 1 matrix where 1 is the class selected and 0 the rejected

ToolBox = getGlobalToolBox;
cArtery = [255 22 18] / 255;
cVein = [18 23 255] / 255;
cChoroid = [0 179 0] / 255;

numClasses = size(classes, 1);
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

[numX, numY] = size(image);
image = rescale(image);
image(isnan(image)) = 0;
level = multithresh(image(mask), numClasses - 1);
graphThreshHistogram(image, level, mask, color, name)
level = [-1 level];
quantizedImage = imquantize(image - 2 * ~mask, level);
quantizedRGB_Mask = repmat(quantizedImage, [1 1 3]);
quantizedImageRGB = zeros(numX, numY, 3);

maskArtery = zeros(numX, numY, 'logical');
maskVein = zeros(numX, numY, 'logical');
maskChoroid = zeros(numX, numY, 'logical');

for i = 1:numClasses
    IDX = reshape([i + 1; i + 1; i + 1], 1, 1, 3);

    if classes(i) == 1
        maskArtery = maskArtery + (quantizedImage == i + 1);
        quantizedImageRGB = quantizedImageRGB + reshape(cArtery, 1, 1, 3) .* (quantizedRGB_Mask == IDX);
    elseif classes(i) == -1
        maskVein = maskVein + (quantizedImage == i + 1);
        quantizedImageRGB = quantizedImageRGB + reshape(cVein, 1, 1, 3) .* (quantizedRGB_Mask == IDX);
    elseif classes(i) == 2
        maskChoroid = maskChoroid + (quantizedImage == i + 1);
        quantizedImageRGB = quantizedImageRGB + reshape(cChoroid, 1, 1, 3) .* (quantizedRGB_Mask == IDX);
    elseif classes(i) == 0
        quantizedImageRGB = quantizedImageRGB + reshape([1 1 1], 1, 1, 3) .* (quantizedRGB_Mask == IDX);
    end

end

imwrite(rescale(quantizedImageRGB), fullfile(ToolBox.path_png, 'mask', 'steps', sprintf("%s_%s_Quantize.png", ToolBox.main_foldername, name)))

maskArtery = logical(maskArtery);
maskVein = logical(maskVein);
maskChoroid = logical(maskChoroid);

maskBackground = ~(mask);
end
