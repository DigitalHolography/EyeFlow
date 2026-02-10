function sectionMontage(subImg_cell, ~, name)

ToolBox = getGlobalToolBox;

numCircles = size(subImg_cell, 1);
numBranches = size(subImg_cell, 2);

rowHeights = zeros(1, numBranches);
colWidths = zeros(1, numCircles);

for i = 1:numCircles
    for j = 1:numBranches
        if ~isempty(subImg_cell{i, j})
            sz = size(subImg_cell{i, j});
            rowHeights(j) = max(rowHeights(j), sz(1));
            colWidths(i) = max(colWidths(i), sz(2));
        end
    end
end

if all(rowHeights == 0) || all(colWidths == 0)
    return
end

rowOffsets = [0 cumsum(rowHeights)];
colOffsets = [0 cumsum(colWidths)];

sub_images_mat = zeros(rowOffsets(end), colOffsets(end));

for i = 1:numCircles
    for j = 1:numBranches
        if ~isempty(subImg_cell{i, j})
            img = subImg_cell{i, j};
            h = size(img, 1);
            w = size(img, 2);
            sub_images_mat( ...
                rowOffsets(j) + (1:h), ...
                colOffsets(i) + (1:w)) = img;
        end
    end
end

sub_images_mat(isnan(sub_images_mat)) = 0;

figure("Visible", "off")
imwrite(sub_images_mat, ...
    fullfile(ToolBox.path_png, ...
    sprintf("%s_%s", ToolBox.folder_name, ...
    sprintf('%s_section_montage.png', name))))

end
