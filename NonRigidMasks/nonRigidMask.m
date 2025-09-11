function nonRigidMask(source, target, aux, warpedAuxPath)
    % Load data
    I1 = imread(source); % reference image (fixed)
    I2 = imread(target); % target image (moving)
    M1 = imread(aux); % mask aligned with I1

    % Convert to grayscale doubles for registration
    I1g = im2double(im2gray(I1));
    I2g = im2double(im2gray(I2));

    % Apply diffeomorphic demons
    [D, M1_warp] = diffeomorphicDemon(I1g, I2g, M1);

    % --- Visualization similar to NonRigidMask ---
    figure("Name", "Rigid+Demons Comparison", "Color", "w");
    tiledlayout(2, 3, "Padding", "compact", "TileSpacing", "compact");

    % Tile 1: Reference
    nexttile; imshow(I1g, []); title("Reference Image");

    % Tile 2: Target
    nexttile; imshow(I2g, []); title("Target Image");

    % Tile 3: Deformation grid visualization
    gridSpacing = 50;
    grid = checkerboard(gridSpacing, ...
        ceil(size(I2g, 1) / (2 * gridSpacing)), ...
        ceil(size(I2g, 2) / (2 * gridSpacing)));
    grid = imresize(grid, size(I2g));
    grid = im2double(grid);
    gridWarped = imwarp(grid, D * 2, "linear");
    nexttile; imshow(gridWarped, []); title("Deformation Grid");

    % Tile 4: Source mask overlay
    nexttile;
    imshow(I1g, []); hold on;
    imshow(labeloverlay(I1g, M1 > 0, 'Transparency', 0.5));
    title("Source Mask on Reference");

    % Tile 5: Warped mask overlay
    tile5 = nexttile;
    imshow(I2g, []); hold on;
    imshow(labeloverlay(I2g, M1_warp > 0, 'Transparency', 0.5));
    title("Warped Mask on Target");
    imwrite(getframe(tile5).cdata, warpedAuxPath)

    % Tile 6: Direct mask comparison
    nexttile; imshowpair(M1, M1_warp); title("Mask Before vs After");

end
