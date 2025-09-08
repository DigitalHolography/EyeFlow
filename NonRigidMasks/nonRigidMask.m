function nonRigidMask(refImg, targImg, maskIn, maskOut)
    I1 = imread(refImg);
    I2 = imread(targImg);
    M1 = imread(maskIn);

    if size(M1,3) == 3
        M1 = rgb2gray(M1);
    end

    % convert image format to rgb
    if size(I1,3) == 3
        I1 = rgb2gray(I1);
    end
    if size(I2,3) == 3
        I2 = rgb2gray(I2); 
    end
    
    I1 = im2double(I1);
    I2 = im2double(I2);

    [optimizer, metric] = imregconfig("monomodal");
    optimizer.MaximumIterations = 100; % this value can be tweaked for either faster or better results
    tformRigid = imregtform(I1, I2, "rigid", optimizer, metric);
    RA = imref2d(size(I1));
    RB = imref2d(size(I2));
    I1_rigid = imwarp(I1, RA, tformRigid, "linear", "OutputView", RB);
    M1_rigid = imwarp(M1,   RA, tformRigid, "nearest", "OutputView", RB);

    diffImg = abs(I2 - I1_rigid);
    if exist("imgradient","file")
        G = imgradient(I2);
    else
        [gx, gy] = gradient(I2);
        G = hypot(gx, gy);
    end
    tauDiff = 3/255;
    tauGrad = 0.02;
    freeze = (diffImg <= tauDiff) & (G <= tauGrad);

    iters = [5 3];
    accFieldSmooth = 10.0;

    [D, ~] = imregdemons( ...
        I1, I2, iters, ...
        "AccumulatedFieldSmoothing", accFieldSmooth, ...
        "PyramidLevels",             numel(iters), ...
        "DisplayWaitbar",            false);

    [M1_warp, ~] = imwarp(M1_rigid, D, "nearest");

    if ~isequal(size(M1_warp), size(I2))
        M1_warp = imresize(M1_warp, size(I2), "nearest");
    end

    M_out = M1_warp;
    M_out(freeze) = M1_rigid(freeze);

    if exist("strel","file") && exist("imclose","file")
        M_out = imclose(M_out, strel("disk",1));
    end

    if nargin < 4 || isempty(maskOut)
        maskOut = "warped_mask.png";
    end
    imwrite(logical(M_out), maskOut);

    gridSpacing = 50;
    grid = checkerboard(gridSpacing, ceil(size(I2,1)/(2*gridSpacing)), ceil(size(I2,2)/(2*gridSpacing)));
    grid = imresize(grid, size(I2));
    grid = im2double(grid);
    G_warped = imwarp(grid, D*5, "linear");
    imwrite(G_warped, "deformationVector.png");

    figure("Name","Rigid+Demons with freeze gating","Color","w");
    tiledlayout(2,3,"Padding","compact","TileSpacing","compact");
    nexttile; imshow(I1,[]);
    title("Source Image")
    nexttile; imshow(I2,[]);
    title("Target Image")
    nexttile; imshow(grid,[]);
    title("Grid (For Reference)")

    % tile 4
    ax4 = nexttile;
    imshow(I1,[],"Parent",ax4); hold(ax4,"on");
    imshow(labeloverlay(I1, M1>0, 'Transparency',0.5), "Parent",ax4);
    title("Source Mask on Source Image")

    % tile 5
    ax5 = nexttile;
    imshow(I2,[],"Parent",ax5); hold(ax5,"on");
    imshow(labeloverlay(I2, M1_warp>0, 'Transparency',0.5), "Parent",ax5);
    title("Mask after Demons on Target Image")

    nexttile; imshow(G_warped,[]);  title("Deformation Vector")
 end