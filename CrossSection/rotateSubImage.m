function [rotatedImg, orientation] = rotateSubImage(subImg, subImgCropped, loc, xy_barycenter)
% Rotate the sub-image to align the blood vessel vertically.
% The orientation is determined by maximizing the kurtosis of the horizontal projection.
%
% Input:
%   subImg - 2D array, the sub-image containing the blood vessel.
%
% Output:
%   rotatedImg - 2D array, the rotated image.
%   orientation - Scalar, the orientation angle used for rotation.

% Initialize an array to store kurtosis values for each angle
subImgCropped(isnan(subImgCropped)) = 0;
subImgCropped(subImgCropped < 0) = 0;
[numX, ~] = size(subImgCropped);
rangeX = floor(numX / 3):ceil(2 * numX / 3);

% Calculate displacement vector components
dx = loc(1) - xy_barycenter(1); % Δx
dy = loc(2) - xy_barycenter(2); % Δy

% Compute angle using atan2 (automatically handles quadrants)
alpha_rad = atan2(dy, dx); % Angle in radians (-π to π)

% Convert to degrees (0° to 360°)
alpha_deg = mod(alpha_rad * 180 / pi, 360);
beta_deg = mod(90 + round(alpha_deg), 180);

% Define a range of angles to test (0 to 180 degrees in 1-degree increments)
N = 40;
angles_deg = (beta_deg - N:beta_deg + N);
angles_deg = mod(angles_deg, 180);
peak_ratio = zeros(1, length(angles_deg));

% Loop over each angle and compute the kurtosis of the horizontal projection
for theta = 1:length(angles_deg)
    % Rotate the image by the current angle
    tmpImg = imrotate(subImgCropped, angles_deg(theta), 'nearest', 'crop');

    % Compute the horizontal projection (sum along rows)
    projx = sum(tmpImg(rangeX, :), 1);

    % Compute kurtosis of the projection (higher kurtosis = more peaked)
    peak_ratio(theta) = max(projx) / sum(projx); % Higher ratio = more spike-like
end

% Find the angle with the maximum kurtosis (most peaked distribution)
[~, idc] = max(peak_ratio);

% Get the corresponding orientation angle
orientation = angles_deg(idc(1));

% Rotate the original image by the determined orientation angle
rotatedImg = imrotatecustom(subImg, orientation);

end
