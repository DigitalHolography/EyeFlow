function saveCrossSectionFigure(subImg, c1, c2, ToolBox, name_section)
% Save a figure showing the cross-section of the blood vessel.
%
% Inputs:
%   subImg              - 2D array, the sub-image of the blood vessel.
%   tmp_section         - 1D array, the cross-section profile.
%   crossSectionWidth   - Scalar, the width of the cross-section.
%   sectionIdx          - Scalar, index of the current section.
%   ToolBox             - Struct, contains parameters and paths.
%   insert              - String, additional identifier for the filename.
%   name_section        - String, name of the section.

% Create figure
f = figure('Visible', 'off');

% Compute Section Cut
[numX, numY] = size(subImg);
profile = sum(subImg, 1, 'omitnan') / numX;
cross_section_profile = (profile ./ max(profile)) * numY;
cross_section_profile(cross_section_profile < 0) = 0;

% Define x-axis values
xAx = linspace(0, numX, numY);

% Display the sub-image
imagesc(xAx, xAx, subImg);
colormap('gray')
axis image;
hold on;

% Flip the y-axis direction
set(gca, 'YDir', 'normal'); % This makes the y-axis increase from bottom to top

% Set aspect ratio
set(gca, 'PlotBoxAspectRatio', [1, 1, 1]);

% Plot the cross-section profile
p = plot(xAx, cross_section_profile);
p.LineWidth = 3;
p.Color = 'red';
p.LineStyle = ':';

% Plot the cross-section width line
x = [c1, c2];
y = [round(numY / 2), round(numY / 2)];
line(x, y, 'Color', 'red', 'LineWidth', 3);
xline(c1, 'r--', 'LineWidth', 3);
xline(c2, 'r--', 'LineWidth', 3);

% Turn off axes
axis off;

% Capture the figure
frame = getframe(gca);

% Save the figure
outputPath = fullfile(ToolBox.path_png, 'vesselSegmentLumen', ...
    sprintf('%s_cross_section_%s.png', ToolBox.folder_name, name_section));
imwrite(frame.cdata, outputPath);

% Close the figure
close(f);
end
