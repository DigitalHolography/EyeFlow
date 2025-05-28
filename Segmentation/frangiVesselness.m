function [M0_binary_img, M0_vesselness_img] = frangiVesselness(M0_ff_img, params)
% Input:
% M0_ff_img - Vessel image (numX x numY)
% params - Structure containing parameters for Frangi filter
%   params.BlackWhite - Use black and white for Frangi filter (default: false)
%   params.range - Frangi scale range (default: [1, 10])
%   params.step - Frangi scale step (default: 1)
% Output:
%   M0_binary_img - Binary mask of segmented vessels
%   M0_vesselness_img - Enhanced vessel image

arguments
    M0_ff_img % Input image (2D grayscale or 3D RGB)
    params.BlackWhite = false % Use black and white for Frangi filter
    params.range = [1, 10] % Frangi scale range
    params.step = 1 % Frangi scale step
end

% Step 1: Preprocessing
% Convert to grayscale if the image is RGB
if size(M0_ff_img, 3) == 3
    M0_gray_img = rgb2gray(M0_ff_img);
else
    M0_gray_img = M0_ff_img;
end

% Normalize the image to [0, 1]
M0_norm_img = rescale(M0_gray_img);

% Step 3: Vessel Enhancement
% Apply Frangi Vesselness Filter
M0_vesselness_img = FrangiFilter2D(M0_norm_img, ...
    "FrangiScaleRange", params.range, ...
    "FrangiScaleRatio", params.step, ...
    "BlackWhite", params.BlackWhite);

% Step 4: Binarization
% Use adaptive thresholding to create a binary mask
M0_vesselness_img(isnan(M0_vesselness_img)) = 0;
M0_binary_img = imbinarize(rescale(M0_vesselness_img), "adaptive");

end
