function [vesselMask, staticVesselMask, filtered_video] = matchedFilterVesselDetection(video, params)
% matchedFilterVesselDetection - Apply matched filter for vessel detection in retinal images
% INPUTS:
%   video: numX x numY x N array of retinal images (grayscale, normalized to [0,1])
%   params: Structure containing parameters for the matched filter
%     .sigma: Standard deviation of Gaussian kernel (default: 1)
%     .length: Length of the kernel (default: 21, should be odd)
%     .angles: Array of angles (in degrees) for rotated kernels (default: 0:15:165)
%     .threshold: Threshold to binarize the filtered response (default: 0.5)
% OUTPUTS:
%   vesselMask: numX x numY x N binary mask indicating vessel presence
%   staticVesselMask: numX x numY binary mask indicating average vessel presence
%   filtered_video: numX x numY x N array of filtered responses
% This function applies a matched filter for vessel detection in a stack of retinal images.

arguments
    video (:, :, :) double {mustBeReal, mustBeNonempty} % 3D array of images
    params.sigma (1, 1) double {mustBePositive} = 1 % Standard deviation of Gaussian kernel
    params.length (1, 1) double {mustBePositive} = 21 % Length of the kernel (should be odd)
    params.angles (1, :) double {mustBeReal, mustBeNonempty} = 0:15:165 % Angles for rotated kernels
    params.threshold (1, 1) double {mustBeReal, mustBeNonnegative} = 0.5 % Threshold for binarization
end

[numX, numY, numFrames] = size(video);
vesselMask = false(numX, numY, numFrames);

% Generate matched filter kernel (Gaussian + negative offset)
kernel = createMatchedFilterKernel(params.sigma, params.length);

video = rescale(video);
filtered_video = zeros(numX, numY, numFrames);

% Process each image in the stack
parfor frameIdx = 1:numFrames
    img = video(:, :, frameIdx);

    maxResponse = zeros(size(img));

    % Apply kernel at multiple orientations
    for angle = params.angles
        rotatedKernel = imrotate(kernel, angle, 'bilinear', 'crop');
        filtered = imfilter(img, rotatedKernel, 'conv', 'same', 'replicate');
        maxResponse = max(maxResponse, filtered);
    end

    % Binarize the response
    vesselMask(:, :, frameIdx) = maxResponse > params.threshold;
    filtered_video(:, :, frameIdx) = maxResponse;
end

staticVesselMask = mean(vesselMask, 3) > params.threshold;

end

% Helper function to create the matched filter kernel
function kernel = createMatchedFilterKernel(sigma, length)
% Create a 1D Gaussian profile
x = linspace(-floor(length / 2), floor(length / 2), length);
gauss1D = exp(-x .^ 2 / (2 * sigma ^ 2));

% Create 2D kernel (Gaussian with negative mean)
kernel = gauss1D' * gauss1D;
kernel = kernel - mean(kernel(:)); % Zero mean (removes DC response)
end
