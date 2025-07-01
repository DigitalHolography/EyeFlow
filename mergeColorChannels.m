function rgbImage = mergeColorChannels(channels, colorMap)
% mergeColorChannels - Merges multiple grayscale channels into a composite RGB image.
%
% Syntax:
%   rgbImage = mergeColorChannels(channels)
%   rgbImage = mergeColorChannels(channels, colorMap)
%
% Inputs:
%   channels  - Cell array of 2D grayscale images (e.g., {ch1, ch2, ch3, ...})
%   colorMap  - (Optional) Nx3 matrix of RGB values (0–1) mapping each channel
%               to a specific color (default: uses predefined color set).
%
% Output:
%   rgbImage  - Output RGB image (same height/width as input channels)

% Validate input
numChannels = numel(channels);
[H, W] = size(channels{1});
for i = 1:3
    channels{i}=squeeze(channels{i});
end
for i = 2:numChannels
    if ~isequal(size(channels{i}), [H, W])
        error('All input channels must have the same dimensions.');
    end
end

% Default color map (up to 7 predefined)
defaultColors = [
    1 0 0;   % Red
    0 1 0;   % Green
    0 0 1;   % Blue
    1 1 0;   % Yellow
    1 0 1;   % Magenta
    0 1 1;   % Cyan
    1 1 1;   % White
    lines;   % in case of high number of channels complet with default lines colormap
    ];

if nargin < 2
    colorMap = defaultColors(1:min(numChannels, size(defaultColors, 1)), :);
elseif size(colorMap, 1) ~= numChannels || size(colorMap, 2) ~= 3
    error('colorMap must be an Nx3 matrix where N matches the number of channels.');
end

% Initialize RGB image
rgbImage = zeros(H, W, 3);

% Merge each channel with color blending
for i = 1:numChannels
    ch = double(channels{i});
    ch = ch / max(ch(:)); % Normalize each channel independently
    for c = 1:3
        rgbImage(:,:,c) = rgbImage(:,:,c) + ch * colorMap(i,c);
    end
end

% Clip values to [0,1]
rgbImage = min(rgbImage, 1);
end
