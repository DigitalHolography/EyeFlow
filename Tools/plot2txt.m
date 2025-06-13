function [] = plot2txt(tabx, taby, filename)
ToolBox = getGlobalToolBox;
tabx = reshape(tabx, [length(tabx) 1]);
taby = reshape(taby, [length(tabx) 1]);
tmp = [tabx, taby];
fileID = fopen(fullfile(ToolBox.path_txt, ...
    sprintf('%s_%s.txt', ToolBox.folder_name, filename)), 'w');
fprintf(fileID, '%f %f \r\n', tmp');
fclose(fileID);
end
