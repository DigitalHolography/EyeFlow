function sectionImage(M0_ff_img, maskLabel, name)

ToolBox = getGlobalToolBox;
[numX, numY] = size(M0_ff_img);
[numCircles, numBranches] = size(maskLabel);

colors = lines(numCircles);
imgRGB = repmat(M0_ff_img, 1, 1, 3);

for cIdx = 1:numCircles
    mask = zeros(numX, numY);

    for bIdx = 1:numBranches
        mask = mask + maskLabel{cIdx, bIdx};
    end

    indxs = find(mask > 0);
    imgRGB(indxs) = colors(cIdx, 1);
    imgRGB(numY * numX + indxs) = colors(cIdx, 2);
    imgRGB(2 * numY * numX + indxs) = colors(cIdx, 3);

    if cIdx > 1 % intersections should be drawn in white
        previous_mask = zeros(numX, numY);

        for bIdx = 1:numBranches
            previous_mask = previous_mask + maskLabel{cIdx - 1, bIdx};
        end

        indxs = find(mask > 0 & previous_mask > 0);
        imgRGB(indxs) = 1;
        imgRGB(numY * numX + indxs) = 1;
        imgRGB(2 * numY * numX + indxs) = 1;
    end

end

figure("Visible", "off")
imshow(imgRGB)
exportgraphics(gca, fullfile(ToolBox.path_png, sprintf("%s_%s.png", ToolBox.folder_name, sprintf('%s_sections%s', name))))

end
