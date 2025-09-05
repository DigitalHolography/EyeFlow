function warp_mask_with_demons(frame1Path, frame2Path, mask1Path, outMaskPath)

    I1 = imread(frame1Path);
    I2 = imread(frame2Path);
    M1 = imread(mask1Path);

    if size(I1,3) == 3, I1 = rgb2gray(I1); end
    if size(I2,3) == 3, I2 = rgb2gray(I2); end
    I1 = im2double(I1);
    I2 = im2double(I2);

    if exist('imhistmatch','file')
        I1m = imhistmatch(I1, I2);
    else
        I1m = I1;
    end

    [optimizer, metric] = imregconfig('monomodal');
    optimizer.MaximumIterations = 100;
    tformRigid = imregtform(I1m, I2, 'rigid', optimizer, metric);
    RA = imref2d(size(I1m));
    RB = imref2d(size(I2));
    I1_rigid = imwarp(I1m, RA, tformRigid, 'linear', 'OutputView', RB);
    M1_rigid = imwarp(M1,   RA, tformRigid, 'nearest', 'OutputView', RB);

    diffImg = abs(I2 - I1_rigid);
    if exist('imgradient','file')
        G = imgradient(I2);
    else
        [gx, gy] = gradient(I2);
        G = hypot(gx, gy);
    end
    tauDiff = 3/255;
    tauGrad = 0.02;
    freeze = (diffImg <= tauDiff) & (G <= tauGrad);

    iters = [60 30 10];
    accFieldSmooth = 2.0;

    [D, movingReg] = imregdemons( ...
        I1, I2, iters, ...
        'AccumulatedFieldSmoothing', accFieldSmooth, ...
        'PyramidLevels',             numel(iters), ...
        'DisplayWaitbar',            false);

    [M1_warp, ~] = imwarp(M1_rigid, D, 'nearest');

    if ~isequal(size(M1_warp), size(I2))
        M1_warp = imresize(M1_warp, size(I2), 'nearest');
    end

    M_out = M1_warp;
    M_out(freeze) = M1_rigid(freeze);

    if exist('strel','file') && exist('imclose','file')
        M_out = imclose(M_out, strel('disk',1));
    end

    if nargin < 4 || isempty(outMaskPath)
        outMaskPath = 'warped_mask.png';
    end
    imwrite(logical(M_out), outMaskPath);



    G = imread("grid.jpg");
    if size(G,3) == 1
        G_resized = imresize(G, size(I2));              % grayscale grid
        G_resized = im2double(G_resized);
        G_warped  = imwarp(G_resized, D, 'linear');     % same size as I2
    else
        % RGB grid: resize each channel and warp as a 3-channel image
        G_resized = im2double(imresize(G, [size(I2,1) size(I2,2)]));
        G_warped  = imwarp(G_resized, D, 'linear');
    end
    imwrite(G_warped, "deformationVector.png");



    figure('Name','Rigid+Demons with freeze gating','Color','w');
    tiledlayout(2,4,'Padding','compact','TileSpacing','compact');
    nexttile; imshow(I1,[]);                   title('Frame 1');
    nexttile; imshow(I2,[]);                   title('Frame 2 (fixed)');
    nexttile; imshowpair(I2, I1_rigid);        title('Rigid align (I1→I2)');
    nexttile; imshow(M1_rigid);                title('Mask @ Frame 1 (rigid)');
    nexttile; imshow(M1_warp);                 title('Mask after Demons');
    nexttile; imshow(M_out);                   title('Final (frozen unchanged)');
    nexttile; imshow(G_warped,[]);  title('Deformation Vector');

    fprintf('Saved warped mask to: %s\n', outMaskPath);
end