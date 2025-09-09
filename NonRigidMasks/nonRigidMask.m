function nonRigidMask(I1, I2, M1, maskOutPath)

    if size(M1,3) == 3
        M1 = rgb2gray(M1);
    end

    % convert image format to gray
    if size(I1,3) == 3
        I1 = rgb2gray(I1);
    end
    if size(I2,3) == 3
        I2 = rgb2gray(I2); 
    end
    
    I1 = im2double(I1);
    I2 = im2double(I2);

    tic;
    [optimizer, metric] = imregconfig("monomodal");
    optimizer.MaximumIterations = 100; % this value can be tweaked for either faster or better results
    tformRigid = imregtform(I1, I2, "rigid", optimizer, metric);
    RA = imref2d(size(I1));
    RB = imref2d(size(I2));
    I1_rigid = imwarp(I1, RA, tformRigid, "linear", "OutputView", RB);
    M1_rigid = imwarp(M1,   RA, tformRigid, "nearest", "OutputView", RB);
    t1 = toc;
    disp(t1)

    tic;
    diffImg = abs(I2 - I1_rigid);
    if exist("imgradient","file")
        G = imgradient(I2);
    else
        [gx, gy] = gradient(I2);
        G = hypot(gx, gy);
    end
    tauDiff = 10/255;
    tauGrad = 0.02;
    freeze = (diffImg <= tauDiff) & (G <= tauGrad);

    iters = [5 15 30];
    accFieldSmooth = 5.0;

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
    t2 = toc;
    disp(t2)

    if exist("strel","file") && exist("imclose","file")
        M_out = imclose(M_out, strel("disk",1));
    end

    if nargin < 4 || isempty(maskOutPath)
        maskOutPath = "warped_mask.png";
    end
    imwrite(logical(M_out), maskOutPath);

    gridSpacing = 50;
    grid = checkerboard(gridSpacing, ceil(size(I2,1)/(2*gridSpacing)), ceil(size(I2,2)/(2*gridSpacing)));
    grid = imresize(grid, size(I2));
    grid = im2double(grid);
    G_warped = imwarp(grid, D, "linear");

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

%{
benchmark results:
for 20 MaximumIterations: gradiant descent 2.9335s, non-rigid updates 0.2519s

for 50 MaximumIterations: gradiant descent 6.2712s, non-rigid updates 0.2471s

for 100 MaximumIterations: gradiant descent 12.1143s, non-rigid updates 0.2593s

for 200 MaximumIterations: gradiant descent 24.3504s, non-rigid updates 0.2485s
%} 