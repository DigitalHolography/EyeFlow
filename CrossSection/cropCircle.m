function img = cropCircle(subImg)
% cropCircle - crops a circular region from a grayscale image using diskMask
%   Input:
%       subImg - Grayscale input image (2D matrix)
%   Output:
%       img - Image with everything outside the circular region set to 0.
%             The circle is centered and has diameter equal to min(numX,numY).

% Get image dimensions
[numX, numY] = size(subImg);

% Create circular mask using diskMask
[mask, ~] = diskMask(numX, numY, 0.5);

% Apply mask to image
img = subImg .* mask;
img(~mask) = NaN;
end
