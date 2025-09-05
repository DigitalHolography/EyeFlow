function warp_mask_with_demons(frame1Path, frame2Path, mask1Path, outMaskPath)

    I1 = imread(frame1Path);
    I2 = imread(frame2Path);
    M1 = imread(mask1Path);

    % convert image format to rgb
    if size(I1,3) == 3
        I1 = rgb2gray(I1);
    end
    if size(I2,3) == 3
        I2 = rgb2gray(I2); 
    end
    
    I1 = im2double(I1);
    I2 = im2double(I2);

    [optimizer, metric] = imregconfig('monomodal');
    optimizer.MaximumIterations = 100;
    tformRigid = imregtform(I1, I2, 'rigid', optimizer, metric);
    RA = imref2d(size(I1));
    RB = imref2d(size(I2));
    I1_rigid = imwarp(I1, RA, tformRigid, 'linear', 'OutputView', RB);
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



    grid = imread("grid.jpg");
    if size(grid,3) == 1
        G_resized = imresize(grid, size(I2));
        G_resized = im2double(G_resized);
        G_warped  = imwarp(G_resized, D, 'linear');
    else
        G_resized = im2double(imresize(grid, [size(I2,1) size(I2,2)]));
        G_warped  = imwarp(G_resized, D, 'linear');
    end
    imwrite(G_warped, "deformationVector.png");



    figure('Name','Rigid+Demons with freeze gating','Color','w');
    tiledlayout(2,3,'Padding','compact','TileSpacing','compact');
    nexttile; imshow(I1,[]);                   title('Frame 1');
    nexttile; imshow(I2,[]);                   title('Frame 2');
    nexttile; imshow(grid,[]);                 title('Grid');
    
    nexttile; imshow(M1);                      title('Mask');
    nexttile; imshow(M1_warp);                 title('Mask after Demons');
    nexttile; imshow(G_warped,[]);  title('Deformation Vector');

    fprintf('Saved warped mask to: %s\n', outMaskPath);
end
