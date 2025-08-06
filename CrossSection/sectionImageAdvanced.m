function sectionImageAdvanced(M0_ff_img, maskLabelArtery, maskLabelVein, maskRejectedArtery, maskRejectedVein, maskVessel)

ToolBox = getGlobalToolBox;
lght_r = [1, 205/255, 210/255]; % light red
red_ = [229/255, 115/255, 115/255]; % red
lght_blu = [179, 229, 252] / 255; % light blue
blu_ = [79, 195, 247] / 255; % blue
[numX, numY] = size(M0_ff_img);
[numCirclesA, numBranchesA] = size(maskLabelArtery);
[numCirclesV, numBranchesV] = size(maskLabelVein);

imgRGB = repmat(M0_ff_img, 1, 1, 3);

% set white the full vessels
indxs = find(maskVessel > 0);
imgRGB(indxs) = 1;
imgRGB(numY * numX + indxs) = 1;
imgRGB(2 * numY * numX + indxs) = 1;

for cIdx = 1:numCirclesA

    if mod(cIdx, 2) == 1
        color = lght_r;
    else
        color = red_;
    end

    mask = zeros(numX, numY);

    for bIdx = 1:numBranchesA
        mask = mask + maskLabelArtery{cIdx, bIdx};
    end

    indxs = find(mask > 0);
    imgRGB(indxs) = color(1);
    imgRGB(numY * numX + indxs) = color(2);
    imgRGB(2 * numY * numX + indxs) = color(3);

    if cIdx > 1 % intersections should be drawn in white
        previous_mask = zeros(numX, numY);

        for bIdx = 1:numBranchesA
            previous_mask = previous_mask + maskLabelArtery{cIdx - 1, bIdx};
        end

        indxs = find(mask > 0 & previous_mask > 0);
        imgRGB(indxs) = 1;
        imgRGB(numY * numX + indxs) = 1;
        imgRGB(2 * numY * numX + indxs) = 1;
    end

end

% note that veins are drawn on top could be fixed
for cIdx = 1:numCirclesV

    if mod(cIdx, 2) == 1
        color = lght_blu;
    else
        color = blu_;
    end

    mask = zeros(numX, numY);

    for bIdx = 1:numBranchesV
        mask = mask + maskLabelVein{cIdx, bIdx};
    end

    indxs = find(mask > 0);
    imgRGB(indxs) = color(1);
    imgRGB(numY * numX + indxs) = color(2);
    imgRGB(2 * numY * numX + indxs) = color(3);

    if cIdx > 1 % intersections should be drawn in white
        previous_mask = zeros(numX, numY);

        for bIdx = 1:numBranchesV
            previous_mask = previous_mask + maskLabelVein{cIdx - 1, bIdx};
        end

        indxs = find(mask > 0 & previous_mask > 0);
        imgRGB(indxs) = 1;
        imgRGB(numY * numX + indxs) = 1;
        imgRGB(2 * numY * numX + indxs) = 1;
    end

end

% set black in the rejected regions
if ~isempty(maskRejectedVein)
    indxs = find(maskRejectedArtery(:, :, 1) | maskRejectedVein(:, :, 1) > 0); % the red channel = rejected parts
else
    indxs = find(maskRejectedArtery(:, :, 1) > 0); % the red channel = rejected parts
end

imgRGB(indxs) = 0;
imgRGB(numY * numX + indxs) = 0;
imgRGB(2 * numY * numX + indxs) = 0;

figure('Visible', 'off');
imshow(imgRGB)
exportgraphics(gca, fullfile(ToolBox.path_png, sprintf("%s_%s.png", ToolBox.folder_name, sprintf('sections_advanced'))))

end
