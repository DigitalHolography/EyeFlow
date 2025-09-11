function [D, warpedAux] = diffeomorphicDemon(source, target, aux)
% diffeomorphicDemon applies rigid + diffeomorphic demons registration
% from 'source' to 'target', and warps 'aux' with this transformation.

% args:
% source: fixed image (double, grayscale in [0,1])
% target: moving image (double, grayscale in [0,1])
% aux: image to which apply the transformation

% returns:
% D: displacement field from source to target
% warpedAux: aux warped into target space

    if size(source,3) == 3, source = rgb2gray(source); end
    if size(target,3) == 3, target = rgb2gray(target); end
    source = im2double(source);
    target = im2double(target);

    % gradient descent
    [optimizer, metric] = imregconfig("monomodal");
    optimizer.MaximumIterations = 100;

    tformRigid = imregtform(source, target, "rigid", optimizer, metric);
    Rtarget = imref2d(size(target));

    sourceRigid = imwarp(source, tformRigid, "linear", "OutputView", Rtarget);
    auxRigid    = imwarp(aux,    tformRigid, "nearest", "OutputView", Rtarget);

    % Diffeomorphic demons
    iters = [30 15 5];
    accSmooth = 10.0;

    [D, ~] = imregdemons(sourceRigid, target, iters, ...
                         "AccumulatedFieldSmoothing", accSmooth, ...
                         "PyramidLevels", numel(iters));

    warpedAux = imwarp(auxRigid, D, "nearest");
end
