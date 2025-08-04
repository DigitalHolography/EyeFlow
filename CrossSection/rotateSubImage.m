function [rotatedImg, orientation] = rotateSubImage(subImg, subImgCropped)
% Rotate the sub-image to align the blood vessel vertically.
% The orientation is determined by maximizing the kurtosis of the horizontal projection.
%
% Input:
%   subImg - 2D array, the sub-image containing the blood vessel.
%
% Output:
%   rotatedImg - 2D array, the rotated image.
%   orientation - Scalar, the orientation angle used for rotation.

% Define a range of angles to test (0 to 180 degrees in 1-degree increments)
angles = linspace(0, 180, 181);

% Initialize an array to store kurtosis values for each angle
peak_ratio = zeros(1, length(angles));
subImgCropped(isnan(subImgCropped)) = 0;
[numX, ~] = size(subImgCropped);
rangeX = floor(numX / 3):ceil(2 * numX / 3);

% Loop over each angle and compute the kurtosis of the horizontal projection
for theta = 1:length(angles)
    % Rotate the image by the current angle
    tmpImg = imrotate(subImgCropped, angles(theta), 'nearest', 'crop');

    % Compute the horizontal projection (sum along rows)
    projx = sum(tmpImg(rangeX, :), 1);

    % Compute kurtosis of the projection (higher kurtosis = more peaked)
    peak_ratio(theta) = max(projx) / sum(projx); % Higher ratio = more spike-like
end

% Find the angle with the maximum kurtosis (most peaked distribution)
[~, idc] = max(peak_ratio);

% Get the corresponding orientation angle
orientation = angles(idc(1));

% Rotate the original image by the determined orientation angle
rotatedImg = imrotatecustom(subImg, orientation);

end
