function [found, diameter_x, diameter_y, x_center, y_center] = findPapilla(M0img)
%Returns the diameter of the papilla measured

try
    if ~isfile('Models/OpticDiscDetector.mat')
        url = 'https://huggingface.co/DigitalHolography/EyeFlow_OpticDiscDetector/resolve/main/OpticDiscDetector.mat';
        websave('Models/OpticDiscDetector.mat', url);
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

data = load('Models/OpticDiscDetector.mat');

M0img = imresize(M0img, [512, 512]);

if size(M0img,3) == 1
    M0img = cat(3, M0img, M0img, M0img);
end

[bboxes, scores, labels] = detect(data.trainedDetector6, M0img, 'Threshold', 0.5);

if ~isempty(bboxes)
    found = true;


    [~, idx] = max(scores);
    bestBox = bboxes(idx, :);

    newWidths = bestBox(:,3) * 0.9;
    newHeights = bestBox(:,4) * 0.9;


    dx = (bestBox(:,3) - newWidths) / 2;
    dy = (bestBox(:,4) - newHeights) / 2;

    newX = bestBox(:,1) + dx;
    newY = bestBox(:,2) + dy;

    shrunkenBBoxes = [newX, newY, newWidths, newHeights];

    diameter_x = shrunkenBBoxes(3);
    diameter_y = shrunkenBBoxes(4);
    x_center = shrunkenBBoxes(1) + diameter_x / 2;
    y_center = shrunkenBBoxes(2) + diameter_y / 2;


    a = shrunkenBBoxes(3)/2;
    b = shrunkenBBoxes(4)/2;

    angle = linspace(0, 2*pi, 100);

    x_ellipsis = x_center + a*cos(angle);
    y_ellipsis = y_center + b*sin(angle);

  
    ToolBox = getGlobalToolBox;

    figure('Visible','off');
    imshow(M0img);
    hold on;
    plot(x_ellipsis, y_ellipsis, 'r', 'LineWidth', 2);
    axis equal;
    exportgraphics(gcf,fullfile(ToolBox.path_png,sprintf('%s_opticdisc.png',ToolBox.folder_name)));
    diameter = (diameter_x + diameter_y) / 2;
    ToolBox.Outputs.add('PapillaRatio',diameter/512,'ua');
    

else
    found = false;
    diameter_x = [94.24];
    diameter_y = [94.24];
    x_center = NaN;
    y_center = NaN;
end


end