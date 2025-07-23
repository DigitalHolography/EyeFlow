function [u, v] = simpleOpticalFlowLK(I1, I2, windowSize)
% SIMPLEOPTICALFLOWLK Estimates optical flow between two images using Lucas-Kanade
% Inputs:
%   I1, I2: Grayscale or binary images (same size)
%   windowSize: Size of local window (odd number, e.g., 5 or 9)
% Outputs:
%   u, v: Horizontal and vertical flow components

% Convert to double
I1 = double(I1);
I2 = double(I2);

% Compute gradients
[Ix, Iy] = gradient(I1);
It = I2 - I1;

% Window (smoothing filter)
w = fspecial('average', windowSize);

% Compute products of derivatives
Ix2 = imfilter(Ix .^ 2, w);
Iy2 = imfilter(Iy .^ 2, w);
Ixy = imfilter(Ix .* Iy, w);
Ixt = imfilter(Ix .* It, w);
Iyt = imfilter(Iy .* It, w);

% Initialize flow
[u, v] = deal(zeros(size(I1)));

% Compute flow (avoid division by zero)
det = (Ix2 .* Iy2 - Ixy .^ 2) +1e-5;
u = (-Iy2 .* Ixt + Ixy .* Iyt) ./ det;
v = (Ixy .* Ixt - Ix2 .* Iyt) ./ det;

end
