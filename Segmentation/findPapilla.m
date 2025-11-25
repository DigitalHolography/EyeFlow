function [found, diameter_x, diameter_y, x_center, y_center] = findPapilla(M0img, net)
%Returns the diameter of the papilla measured
ToolBox = getGlobalToolBox;
params = ToolBox.getParams;
saveFigures = params.saveFigures;

% Initialize output variables
[numX, numY] = size(M0img);
found = false;
diameter_x = NaN;
diameter_y = NaN;
x_center = numX / 2;
y_center = numY / 2;

if ~params.json.Mask.OpticDiskDetectorNet
    return;
end

if size(M0img, 3) == 1
    M0img = repmat(M0img, 1, 1, 3);
end

M0img = imresize(rescale(M0img), [512, 512]); % 512 is the iniput size of the net
input = rescale(M0img);
output = predict(net, input);

if ~isempty(output)
    % Find the box with the highest confidence
    [val, idx] = max(output(:, 5, :));
    bestBox = output(:, :, idx);

    % Check if the confidence is above a threshold (e.g., 0.5)
    if val > 0.5
        found = true;
        diameter_x = bestBox(3);
        diameter_y = bestBox(4);
        x_center = round(bestBox(1) / 512 * numX);
        y_center = round(bestBox(2) / 512 * numY);
        ToolBox.Output.DimOut.add("PapillaRatio", (diameter_x + diameter_y) / 2/512, ["test"], '');
    end

    if saveFigures
        a = bestBox(3) / 2;
        b = bestBox(4) / 2;
        angle = linspace(0, 2 * pi, 360);
        x_ellipsis = x_center * 512 / numX + a * cos(angle);
        y_ellipsis = y_center * 512 / numY + b * sin(angle);

        figure('Visible', 'off', 'Color', 'w');
        imshow(M0img, []);
        hold on;

        if val > 0.5
            plot(x_ellipsis, y_ellipsis, 'g--', 'LineWidth', 2);
        else
            plot(x_ellipsis, y_ellipsis, 'r--', 'LineWidth', 2);
        end

        axis equal;
        axis off;
        exportgraphics(gcf, fullfile(ToolBox.path_png, sprintf('%s_opticdisc.png', ToolBox.folder_name)));
    end

end

close all

end
