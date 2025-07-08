function [warped , ux, uy] = blocWiseAffineRegistration(moving,fixed,NumPatch,SmoothingValue)
% bloc by bloc displacement field estimation using only translational local
% displacement estimate (registerImagesCrossCorrelation) and gaussian smoothing
%  

[H, W] = size(fixed);
ux = zeros(H, W);
uy = zeros(H, W);
patchSize = floor(max(H,W)/NumPatch);

% Init affine optimizer
% [optimizer, metric] = imregconfig('multimodal');
% metric = registration.metric.MeanSquares;

% Local affine registration per patch
for y = 1:patchSize:H - patchSize
    for x = 1:patchSize:W - patchSize
        % Extract patches
        fixedPatch = fixed(y:y+patchSize-1, x:x+patchSize-1);
        movingPatch = moving(y:y+patchSize-1, x:x+patchSize-1);

        % Get coordinates of patch
        % [Xp, Yp] = meshgrid(x:x+patchSize-1, y:y+patchSize-1);
        % coords = [Xp(:), Yp(:)];

        % Estimate affine transform (e.g., using intensity-based optimization or features)
        %tform = imregtform(rescale(movingPatch), rescale(fixedPatch), ...
        %'affine',optimizer, metric, 'DisplayOptimization', false); nul
        [~, shift] = registerImagesCrossCorrelation(movingPatch, fixedPatch);
        reg_image = circshift(movingPatch, shift);
        % shift
        % Transform coordinates
        if ~isempty(shift) %~isempty(tform)
            % newCoords = transformPointsForward(tform, coords);
            dx = (shift(2)+patchSize)*ones(patchSize);
            dy = (shift(1)+patchSize)*ones(patchSize);

            % Accumulate displacements
            ux(y:y+patchSize-1, x:x+patchSize-1) = dx;
            uy(y:y+patchSize-1, x:x+patchSize-1) = dy;
        end
        % figure,imshowpair(fixedPatch,reg_image);

    end
end

% Smooth the displacement field

if SmoothingValue > 0 % usually 2
    ux = imgaussfilt(ux, SmoothingValue);
    uy = imgaussfilt(uy, SmoothingValue);
end

% Warp the image
[Xq, Yq] = meshgrid(1:W, 1:H);
warped = interp2(moving, Xq - ux, Yq - uy, 'linear', 0);

imshowpair(fixed, warped);

end