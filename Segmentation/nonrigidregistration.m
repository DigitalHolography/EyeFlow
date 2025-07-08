function [warpedMask,similarMaskAffine,warpedMaskDilated] = nonrigidregistration(similarM0,M0,similarMask,folderpath,tag)
% Non rigid registration based on non rigid algotrithms (Bloc wise
% registration for now)
% In our case M0 is the current M0 and similarM0 is the old M0 for which we
% have a good segmentation similarMask that we want to project to the
% registration 

fprintf("Wrapping a similar correct mask to the non-rigid registration calculated with non rigid algorithms.\n");

% Step 1: Correct illumination and renormalize within the diaphragm
[numX,numY] = size(similarM0);
disk = diskMask(numX,numY,1/2-0.05);

M0 = rescale(M0);
similarM0 = rescale(similarM0);

M0Diaph = (M0 .*disk - mean(M0(disk)))/max(abs(M0(disk))); % version normalized and with zeros outside the diaphragm
similarM0Diaph = similarM0 .*disk - mean(similarM0(disk))/max(abs(similarM0(disk))); % version normalized and with zeros outside the diaphragm

%M0 = imhistmatch(M0,similarM0);

% Step 1: Compute affine registration 
fprintf("Calculating affine registration : \n");tic;

[optimizer, metric] = imregconfig('multimodal');
tform = imregtform(M0Diaph, similarM0Diaph, 'affine',optimizer, metric, 'DisplayOptimization', false);
tformInv = invert(tform);
Rfixed = imref2d(size(similarM0));
M0Affine = imwarp(M0, tform, 'OutputView', Rfixed);
similarMaskAffine = imwarp(similarMask, tformInv, 'OutputView', Rfixed);
toc;
figure(51),imshowpair(similarM0Diaph,M0Diaph);
saveas(gcf,fullfile(folderpath,sprintf('%s_1_Diff.png',tag)));
figure(52),imshowpair(similarM0,M0Affine);
saveas(gcf,fullfile(folderpath,sprintf('%s_1_AfterAffine.png',tag)));
figure(53),imshowpair(M0,similarMaskAffine);
saveas(gcf,fullfile(folderpath,sprintf('%s_1_MaskAfterAffine.png',tag)));

% Step 2: Compute deformation field using blocwiseAffineRegistration
% 
% "the Demons algorithm
% % % fprintf("Calculating non rigid registration : \n");tic;
% % % [displacementField, ~] = imregdemons(toWarpAffine, similarM0Diaph, ... % Here it could potentially be M0Diaph and similarM0Diaph
% % %     300, ...                    % Number of iterations
% % %     'PyramidLevels', 3, ...
% % %     'AccumulatedFieldSmoothing', 1.3);  % Regularization
% % % 
% % % % Step 3: Warp the source mask using the displacement field
% % % % Create a displacement field in a format suitable for imwarp
% % % D = cat(3, displacementField(:,:,2), displacementField(:,:,1));  % [X,Y] order
% % % 
% % % % Apply deformation using imwarp
% % % % tform = affine2d(eye(3));  % Identity transform
% % % % RA = imref2d(size(M0));  % Reference frame for output
% % % % 'D' is the displacement field (in [X, Y] order)
% % % warpedMask = imwarp(toWarpAffine, +1*D, "nearest");
% warpedMask = toWarpAffine; 
NUM_BLOCS = 32;
SMOOTH_PIX = 1;
M0Frangi = FrangiFilter2D(M0Affine);
similarM0Frangi = FrangiFilter2D(mat2gray(similarM0));
M0FrangiDiaph = (M0Frangi .*disk - mean(M0Frangi(disk)))/max(abs(M0Frangi(disk))); % version normalized and with zeros outside the diaphragm
similarM0FrangiDiaph = (similarM0Frangi .*disk - mean(similarM0Frangi(disk)))/max(abs(similarM0Frangi(disk))); % version normalized and with zeros outside the diaphragm

[warped , u, v] = blocWiseAffineRegistration(M0FrangiDiaph,similarM0FrangiDiaph,NUM_BLOCS,SMOOTH_PIX);
[Xq, Yq] = meshgrid(1:numX, 1:numY);
warpedMask = interp2(single(similarMaskAffine), Xq + u, Yq + v, 'linear', 0)>0;

figure(55),imshowpair(M0,warpedMask);
saveas(gcf,fullfile(folderpath,sprintf('%s_2_MaskBlocWiseReg.png',tag)));
figure(56),imshowpair(FrangiFilter2D(mat2gray(similarM0Diaph)), warped);
saveas(gcf,fullfile(folderpath,sprintf('%s_2_BlocWiseRegWarped.png',tag)));


toc;


step = 50;  % Adjust to make arrows sparser or denser
[xGrid, yGrid] = meshgrid(1:step:size(u,2), 1:step:size(u,1));
uSampled = - u(1:step:end, 1:step:end);
vSampled = - v(1:step:end, 1:step:end);
figure(57); % Plot quiver over mask (optional background)
imshow(similarMaskAffine); hold on;
quiver(xGrid, yGrid, uSampled, vSampled, 'r');  % red arrows
title('Displacement Field (Quiver Plot)');
axis on;
saveas(gcf,fullfile(folderpath,sprintf('%s_2_DisplacementFieldNonRigid.png',tag)));


% Step 4: Dilate the warped mask slightly
SE = strel('disk', 1 ...
    );
warpedMaskDilated = imdilate(warpedMask,SE);

% warpedMaskDilated = imdilate(warpedMaskDilated,SE);
% warpedMaskDilated = imerode(warpedMaskDilated,SE);


figure(58),imshowpair(M0,warpedMaskDilated);
saveas(gcf,fullfile(folderpath,sprintf('%s_3_MaskAfterDilation.png',tag)));


close all;

end