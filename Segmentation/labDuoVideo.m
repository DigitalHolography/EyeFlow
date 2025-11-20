function [RGB_video] = labDuoVideo(video_1, image_2, h)
% labDuoVideo
% Applies labDuoImage-like colorization to each frame of a video
%
% INPUTS:
%   video_1 : 3D or 4D array (HxWxN or HxWx1xN) representing grayscale video frames
%   image_2 : 2D array (HxW) providing color modulation
%   h       : hue angle in degrees (default = 45)
%
% OUTPUT:
%   RGB_video : 4D array (HxWx3xN) containing colorized video frames

arguments
    video_1
    image_2
    h = 45
end

% Normalize hue
h = mod(h, 360);

% Compute rx, ry according to hue
if abs(h) <= 45 || abs(h - 180) <= 45
    rx = sign(cosd(h));
    ry = sign(sind(h)) .* (1 / (cosd(h) * cosd(h)) - 1);
else
    ry = sign(sind(h));
    rx = sign(cosd(h)) .* (1 / (sind(h) * sind(h)) - 1);
end

% Normalize image_2
image_2 = image_2 ./ max(abs(max(image_2, [], "all")), abs(min(image_2, [], "all")));

% Get video dimensions
[H, W, ~, N] = size(video_1);

if ndims(video_1) == 3
    N = size(video_1, 3);
end

% Initialize output
RGB_video = zeros(H, W, 3, N);

% Process each frame
for i = 1:N
    frame = video_1(:, :, i);

    % Normalize frame
    frame = frame ./ max(abs(max(frame, [], "all")), abs(min(frame, [], "all")));

    % Build Lab channels
    L = 100 .* frame;
    a = 128 .* image_2 * rx;
    b = 128 .* image_2 * ry;

    Lab = cat(3, L, a, b);
    RGB_video(:, :, :, i) = lab2rgb(Lab);
end

end
