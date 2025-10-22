function [D, warpedAux] = diffeomorphicDemon(source, target, aux)
% diffeomorphicDemon applies rigid + diffeomorphic demons registration
% from 'source' to 'target', and warps 'aux' with this transformation.

% args:
% source: moving image (double, grayscale in [0,1])
% target: fixed image (double, grayscale in [0,1])
% aux: image to which apply the transformation

% returns:
% D: displacement field from source to target
% warpedAux: aux warped into target
if size(source, 3) == 3, source = rgb2gray(source); end
if size(target, 3) == 3, target = rgb2gray(target); end
source = im2double(source);
target = im2double(target);

% Diffeomorphic demons
iters = [5];
accSmooth = 1.0;

[D, ~] = imregdemons(source, target, iters, ...
    "AccumulatedFieldSmoothing", accSmooth, ...
    "PyramidLevels", numel(iters), ...
    "DisplayWaitbar", false);

warpedAux = imwarp(aux, D, "nearest");

%freeze pixel where warp is minimal
mask = hypot(D(:, :, 1), D(:, :, 2)) < 0.5;
warpedAux(mask) = source(mask);
end
