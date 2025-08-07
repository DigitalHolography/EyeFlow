function [found, diameter_x, diameter_y, x_center, y_center] = findPapilla(M0img)
%Returns the diameter of the papilla measured
ToolBox = getGlobalToolBox;
params = ToolBox.getParams;
mask_params = params.json.Mask;

if mask_params.OpticDiscDetector

    try

        if ~isfile('Models/opticdisc.onnx')
            url = 'https://huggingface.co/DigitalHolography/EyeFlow_OpticDiscDetectorV2/resolve/main/opticdisc.onnx';
            websave('Models/opticdisc.onnx', url);
        end

    catch ME
        disp(ME)
        found = false;
        diameter_x = NaN;
        diameter_y = NaN;
        x_center = NaN;
        y_center = NaN;
        return
    end

    net = importNetworkFromONNX('Models/opticdisc.onnx');

    if size(M0img, 3) == 1
        M0img = repmat(M0img, 1, 1, 3);
    end

    M0img = imresize(rescale(M0img), [512, 512]);

    input = rescale(M0img);

    output = predict(net, input);

    if ~isempty(output)
        found = true;
        [val, idx] = max(output(:, 5, :));
        bestBox = output(:, :, idx);
        % bestScore = val;
        % bestLabel = 'optic disc';

        diameter_x = bestBox(3);
        diameter_y = bestBox(4);

        x_center = bestBox(1);
        y_center = bestBox(2);
        a = bestBox(3) / 2;
        b = bestBox(4) / 2;

        angle = linspace(0, 2 * pi, 360);
        x_ellipsis = x_center + a * cos(angle);
        y_ellipsis = y_center + b * sin(angle);

        ToolBox = getGlobalToolBox;

        figure('Visible', 'off');
        imshow(M0img, []);
        hold on;

        if val > 0.5
            plot(x_ellipsis, y_ellipsis, 'g--', 'LineWidth', 2);
        else
            plot(x_ellipsis, y_ellipsis, 'r--', 'LineWidth', 2);
            found = false;
            diameter_x = [94.24];
            diameter_y = [94.24];
            x_center = NaN;
            y_center = NaN;
        end

        % imgOut = insertObjectAnnotation(M0img, 'rectangle', bestBox, ...
        % sprintf('%s: %.2f', string(bestLabel), bestScore));
        % imshow(imgOut);
        axis equal;
        axis off;
        exportgraphics(gcf, fullfile(ToolBox.path_png, sprintf('%s_opticdisc.png', ToolBox.folder_name)));

    else
        found = false;
        diameter_x = [94.24];
        diameter_y = [94.24];
        x_center = NaN;
        y_center = NaN;
    end

else
    found = false;
    diameter_x = [94.24];
    diameter_y = [94.24];
    x_center = NaN;
    y_center = NaN;
end

end
