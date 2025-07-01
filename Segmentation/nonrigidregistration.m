function [warpedMask,toWarpAffine,warpedMaskDilated] = nonrigidregistration(imFixed,imMoving,prevMask,toWarpMask,folderpath,tag)
% Non rigid registration based on Demons algotrithm
% In our case fixed is the current M0 and moving is the old M0 for which we
% have a good segmentation toWarpMask that we want to project to the
% registration using a preMask prevMask.

fprintf("Wrapping a similar correct mask to the non-rigid registration calculated with Demons algorithm.\n");

% Step 1: Correct illumination and renormalize within the diaphragm
[numX,numY] = size(imFixed);
disk = diskMask(numX,numY,1/2-0.05);

imMoving = rescale(imMoving);
imFixed = rescale(imFixed);

imMovingDiaph = (imMoving .*disk - mean(imMoving(disk)))/max(abs(imMoving(disk))); % version normalized and with zeros outside the diaphragm
imFixedDiaph = imFixed .*disk - mean(imFixed(disk))/max(abs(imFixed(disk))); % version normalized and with zeros outside the diaphragm

%imMoving = imhistmatch(imMoving,imFixed);

% Step 1: Compute affine registration 
fprintf("Calculating affine registration : \n");tic;

[optimizer, metric] = imregconfig('monomodal');
tform = imregtform(imMovingDiaph, imFixedDiaph, 'affine',optimizer, metric);
Rfixed = imref2d(size(imFixed));
movingAffine = imwarp(imMoving, tform, 'OutputView', Rfixed);
toWarpAffine = imwarp(toWarpMask, tform, 'OutputView', Rfixed);
toc;

figure(51),imshowpair(imFixedDiaph,imMovingDiaph);
saveas(gcf,fullfile(folderpath,sprintf('%s_SimilarAffine.png',tag)));
figure(52),imshowpair(imFixed,movingAffine);
saveas(gcf,fullfile(folderpath,sprintf('%s_SimilarAfterAffine.png',tag)));

figure(53),imshowpair(imFixed,toWarpAffine);
saveas(gcf,fullfile(folderpath,sprintf('%s_SimilarMaskAfterAffine.png',tag)));

% Step 2: Compute deformation field using the Demons algorithm
fprintf("Calculating non rigid registration : \n");tic;
[displacementField, ~] = imregdemons(toWarpAffine, prevMask, ... % Here it could potentially be imMovingDiaph and imFixedDiaph
    100, ...                    % Number of iterations
    'AccumulatedFieldSmoothing', 1.3);  % Regularization

% Step 3: Warp the source mask using the displacement field
% Create a displacement field in a format suitable for imwarp
D = cat(3, displacementField(:,:,2), displacementField(:,:,1));  % [X,Y] order

% Apply deformation using imwarp
% tform = affine2d(eye(3));  % Identity transform
% RA = imref2d(size(imMoving));  % Reference frame for output
% 'D' is the displacement field (in [X, Y] order)
warpedMask = imwarp(toWarpAffine, +1*D, "nearest");
toc;

figure(54),imshowpair(imFixed,warpedMask);
saveas(gcf,fullfile(folderpath,sprintf('%s_SimilarMaskAfterDemons.png',tag)));

figure(55),imshowpair(prevMask,warpedMask);
saveas(gcf,fullfile(folderpath,sprintf('%SimilarMaskAfterDemonsOnAuto.png',tag)));


u = displacementField(:,:,2);  % X (horizontal)
v = displacementField(:,:,1);  % Y (vertical)
step = 50;  % Adjust to make arrows sparser or denser
[xGrid, yGrid] = meshgrid(1:step:size(u,2), 1:step:size(u,1));
uSampled = u(1:step:end, 1:step:end);
vSampled = v(1:step:end, 1:step:end);
figure(56); % Plot quiver over mask (optional background)
imshow(toWarpAffine); hold on;
quiver(xGrid, yGrid, uSampled, vSampled, 'r');  % red arrows
title('Displacement Field (Quiver Plot)');
axis on;

% Step 4: Dilate the warped mask slightly on the prevMask
SE = strel('disk', 1);
warpedMaskDilated = warpedMask | prevMask & imdilate(warpedMask,SE);

figure(57),imshowpair(imFixed,warpedMaskDilated);
saveas(gcf,fullfile(folderpath,sprintf('%s_SimilarMaskAfterDilation.png',tag)));

figure(58),imshowpair(prevMask,warpedMaskDilated);
saveas(gcf,fullfile(folderpath,sprintf('%s_SimilarMaskAfterDilationOnAuto.png',tag)));


close all;

end