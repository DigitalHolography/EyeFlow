function blocWiseAffineRegistration(moving,fixed,patchSize)

[H, W] = size(fixed);
ux = zeros(H, W);
uy = zeros(H, W);

% Init affine optimizer
[optimizer, metric] = imregconfig('monomodal');
metric = registration.metric.MeanSquares;

% Local affine registration per patch
for y = 1:patchSize:H - patchSize
    for x = 1:patchSize:W - patchSize
        % Extract patches
        fixedPatch = fixed(y:y+patchSize-1, x:x+patchSize-1);
        movingPatch = moving(y:y+patchSize-1, x:x+patchSize-1);

        % Get coordinates of patch
        [Xp, Yp] = meshgrid(x:x+patchSize-1, y:y+patchSize-1);
        coords = [Xp(:), Yp(:)];

        % Estimate affine transform (e.g., using intensity-based optimization or features)
        tform = imregtform(movingPatch, fixedPatch, 'affine',optimizer, metric, 'DisplayOptimization', false);

        % Transform coordinates
        if ~isempty(tform)
            newCoords = transformPointsForward(tform, coords);
            dx = reshape(newCoords(:,1) - coords(:,1), patchSize, patchSize);
            dy = reshape(newCoords(:,2) - coords(:,2), patchSize, patchSize);

            % Accumulate displacements
            ux(y:y+patchSize-1, x:x+patchSize-1) = dx;
            uy(y:y+patchSize-1, x:x+patchSize-1) = dy;
        end
    end
end

% Smooth the displacement field
% ux = imgaussfilt(ux, 2);
% uy = imgaussfilt(uy, 2);

% Warp the image
[Xq, Yq] = meshgrid(1:W, 1:H);
warped = interp2(moving, Xq + ux, Yq + uy, 'linear', 0);

imshowpair(fixed, warped);

end