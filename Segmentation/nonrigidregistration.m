function [warpedMask] = nonrigidregistration(maskSource,maskTarget)
% Non rigid registration based on Demons algotrithm
% Step 1: Load binary masks (logical arrays)

% Step 2: Compute deformation field using the Demons algorithm
[displacementField, movingReg] = imregdemons(maskSource, maskTarget, ...
    100, ...                    % Number of iterations
    'AccumulatedFieldSmoothing', 1.0);  % Regularization

% Step 3: Warp the source mask using the displacement field
% Create a displacement field in a format suitable for imwarp
D = cat(3, displacementField(:,:,2), displacementField(:,:,1));  % [X,Y] order

% Apply deformation using imwarp
tform = affine2d(eye(3));  % Identity transform
RA = imref2d(size(maskTarget));  % Reference frame for output
warpedMask = imwarp(maskSource, D, tform, 'OutputView', RA, 'InterpolationMethod', 'nearest');

% Step 4: Visualize results
figure;
subplot(1, 3, 1);
imshow(maskSource);
title('Original Mask');

subplot(1, 3, 2);
imshow(maskTarget);
title('Target Mask');

subplot(1, 3, 3);
imshow(warpedMask);
title('Warped Mask (Deformed)');
end