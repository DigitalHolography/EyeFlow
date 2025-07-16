function sectionMontage(subImg_cell, numSections, name)

ToolBox = getGlobalToolBox;

subImgSize = [];

for i = 1:size(subImg_cell, 1)

    for n = 1:size(subImg_cell, 2)

        if ~isempty(subImg_cell{i, n})
            subImgSize = size(subImg_cell{i, n});
        end

        if ~isempty(subImgSize)
            break;
        end

    end

    if ~isempty(subImgSize)
        break;
    end

end

if isempty(subImgSize)
    return
end

figure("Visible", "off")
numCircles = size(subImg_cell, 1);
numBranches = size(subImg_cell, 2);
% fill with zero images the zeros parts

sub_images_mat = zeros((subImgSize(1) * numBranches) + 1, (subImgSize(2) * numCircles) + 1);

for circleIdx = 1:numCircles

    for branchIdx = 1:numBranches

        if ~isempty(subImg_cell{circleIdx, branchIdx})
            sub_images_mat( ...
                1 + (branchIdx - 1) * subImgSize(1):branchIdx * subImgSize(1), ...
                1 + (circleIdx - 1) * subImgSize(2):circleIdx * subImgSize(2)) = subImg_cell{circleIdx, branchIdx};
        end

    end

end

sub_images_mat(isnan(sub_images_mat)) = 0;
imwrite(sub_images_mat, fullfile(ToolBox.path_png, 'local', ...
    sprintf("%s_%s", ToolBox.folder_name, sprintf('%s_section_montage.png', name))))

end
