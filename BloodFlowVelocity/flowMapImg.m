function [combined_img] = flowMapImg(img, masks, cmaps, opt)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

arguments
    img % input image
    masks % cell array of logical masks
    cmaps % cell array of colormaps
    opt.circles = [] % logical mask for circles
    opt.background = [] % background image
end

N = size(masks, 2);

if N ~= size(cmaps, 2)
    error('Number of masks and colormaps must be the same')
end

if ~isempty(opt.circles)

    for i = 1:N
        masks{N + i} = masks{i} & opt.circles;
        masks{i} = masks{i} & ~opt.circles;
        cmaps{N + i} = flip(cmaps{i}, 1);
    end

    N = 2 * N;
end

[numX, numY, ~] = size(img);

layers = zeros(numX, numY, 3, N, 'single');

for i = 1:size(masks, 2)
    mask = masks{i};
    cmap = cmaps{i};
    layers(:, :, :, i) = setcmap(img, mask, cmap);
end

maskVessel = zeros(numX, numY, 'logical');

for i = 1:N
    maskVessel = maskVessel + masks{i};
end

if isempty(opt.circles)
    opt.circles = zeros(numX, numY, 'logical');
end

if isempty(opt.background)
    opt.background = zeros(numX, numY, 'single');
end

combined_img = sum(layers, 4) + (opt.background .* ~maskVessel .* ~opt.circles) + opt.circles .* ~maskVessel;

combined_img = uint8(rescale(combined_img, 0, 255));

end
