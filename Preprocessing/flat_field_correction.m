function corrected_image = flat_field_correction(image, correctionParams, borderAmount, method)
% FLAT_FIELD_CORRECTION Corrects an image for uneven illumination using either Gaussian blur or fitted Gaussian.
%
% Inputs:
%   image: Input image (2D array).
%   correctionParams: Parameters for the correction method.
%       - For 'gaussianBlur': gw (Gaussian blur width).
%       - For 'fittedGaussian': [A, mu, sigma, C] (Gaussian parameters).
%   borderAmount: Fraction of the image border to exclude (optional, default = 0).
%   method: Correction method to use ('gaussianBlur' or 'fittedGaussian').
%
% Output:
%   corrected_image: Flat-field corrected image.

% Set default borderAmount if not provided
if nargin < 3
    borderAmount = 0;
end

% Set default method if not provided
if nargin < 4
    method = 'gaussianBlur'; % Default to Gaussian blur method
end

% Check if image is normalized between 0 and 1
Im_min = min(image, [], 'all');
Im_max = max(image, [], 'all');

if Im_min < 0 || Im_max > 1

    if Im_max > Im_min
        image = (image - Im_min) ./ (Im_max - Im_min);
    else
        image = zeros(size(image), 'like', image); % constant image
    end

    flag = 1;
else
    flag = 0;
end

% Define the non-border region
if borderAmount == 0
    a = 1;
    b = size(image, 1);
    c = 1;
    d = size(image, 2);
else
    a = ceil(size(image, 1) * borderAmount);
    b = floor(size(image, 1) * (1 - borderAmount));
    c = ceil(size(image, 2) * borderAmount);
    d = floor(size(image, 2) * (1 - borderAmount));
end

% Compute the sum of intensities in the non-border region
ms = sum(image(a:b, c:d), [1 2]);

% Apply the selected correction method
switch method
    case 'gaussianBlur'
        % Gaussian blur method
        gw = correctionParams; % Gaussian blur width
        image = image ./ imgaussfilt(image, gw);

    case 'fittedGaussian'
        % Fitted Gaussian method
        gaussianParams = correctionParams; % Gaussian parameters [A, mu, sigma, C]
        [numX, numY, ~] = size(image);
        [X, Y] = meshgrid(1:numY, 1:numX);
        R = sqrt((X - numX / 2) .^ 2 + (Y - numY / 2) .^ 2); % Distance from the center
        illumination_pattern = gaussianParams(1) * exp(- ((R) / gaussianParams(3)) .^ 2) + gaussianParams(4);
        image = image ./ illumination_pattern;

    otherwise
        error('Invalid method. Choose "gaussianBlur" or "fittedGaussian".');
end

% Rescale to maintain total intensity in the non-border region
ms2 = sum(image(a:b, c:d), [1 2]);
corrected_image = (ms / ms2) .* image;

% Rescale back to original intensity range if normalization was applied
if flag == 1
    corrected_image = Im_min + (Im_max - Im_min) .* corrected_image;
end

end
