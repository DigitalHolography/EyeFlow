function sectionMontage(subImg_cell, numSections, name)

ToolBox = getGlobalToolBox;

subImgSize = [];
i = 1;

while isempty(subImgSize)

    if ~isempty(subImg_cell{1, i})
        subImgSize = size(subImg_cell{1, i}{1, 1});
    end

    i = i + 1;
end

numCircles = size(subImg_cell, 2);
figure("Visible", "off")
numSectionMax = max(numSections);
% fill with zero images the zeros parts

sub_images_mat = zeros(subImgSize(1), subImgSize(2), numSectionMax * numCircles);

for circleIdx = 1:numCircles

    for sectionIdx = 1:numSectionMax

        if size(subImg_cell{1, circleIdx}, 2) < numSectionMax
            subImg_cell{1, circleIdx}(end + 1) = {zeros(subImgSize, 'single')};
        end

        sub_images_mat(:, :, ((sectionIdx - 1) * numCircles) + circleIdx) = subImg_cell{1, circleIdx}{1, sectionIdx};

    end

end

montage(sub_images_mat, "Size", [max(1, numSectionMax), max(1, numCircles)]);
exportgraphics(gca, fullfile(ToolBox.path_png, 'crossSectionsAnalysis', sprintf("%s_%s", ToolBox.main_foldername, sprintf('%s_section_montage.png', name))))

end
