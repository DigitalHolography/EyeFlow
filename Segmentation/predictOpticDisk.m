function [opticDiskMask, cx, cy, width, height] = predictOpticDisk(model, M0)
% predictOpticDisk - Predicts the optic disk of a M0 image
%   Inputs:
%       model     - The pre-loaded YOLO model
%       M0        - Input image (2D matrix)
%   Output:
%       opticDiskMask - Logical matrix (mask)
%       cx            - Ellipse center x
%       cy            - Ellipse center y
%       width         - Ellipse width
%       height        - Ellipse height

% Preprocess
inputSize = [1024, 1024]; 
[origH, origW] = size(M0);
resizedM0 = rescale(imresize(M0, inputSize));
rgbM0 = cat(3, resizedM0, resizedM0, resizedM0);
dlInput = dlarray(rgbM0, 'SSCB');

% Predict
[out1, ~] = predict(model, dlInput);

 % --- 3. Extract and Parse Data ---
raw_preds = extractdata(gather(out1));
predictions = squeeze(raw_preds)'; 

% --- 4. Filter and NMS ---
confThreshold = 0.25;
scores = predictions(:, 5);

validIdx = scores > confThreshold;
validPreds = predictions(validIdx, :);

if isempty(validPreds)
    opticDiskMask = false(origH, origW);
    cx = 0; cy = 0; width = 0; height = 0;
    imwrite(rgbM0, 'predicted_optic_disk.png'); % Save blank if no detection
    return;
end

boxes_cxcywh = validPreds(:, 1:4);

% Convert to Top-Left for NMS
boxes_tlwh = boxes_cxcywh;
boxes_tlwh(:, 1) = boxes_cxcywh(:, 1) - boxes_cxcywh(:, 3) / 2;
boxes_tlwh(:, 2) = boxes_cxcywh(:, 2) - boxes_cxcywh(:, 4) / 2;

% Select strongest bounding box
[~, ~, idx] = selectStrongestBbox(boxes_tlwh, validPreds(:, 5), 'OverlapThreshold', 0.6);

bestIdx = idx(1);
bestBox = boxes_cxcywh(bestIdx, :); % [cx, cy, w, h] relative to 1024x1024

% --- 5. Scale Box Coordinates ---
% Scale factors
scaleX = origW / inputSize(2);
scaleY = origH / inputSize(1);

cx = bestBox(1) * scaleX;
cy = bestBox(2) * scaleY;
width = bestBox(3) * scaleX;
height = bestBox(4) * scaleY;

% --- 6. Generate Output Mask (Geometric Ellipse) ---
% Create a coordinate grid for the original image size
[X, Y] = meshgrid(1:origW, 1:origH);

% Ellipse radii
rx = width / 2;
ry = height / 2;

% Standard Ellipse Equation: ((x-h)^2 / rx^2) + ((y-k)^2 / ry^2) <= 1
if rx > 0 && ry > 0
    normDist = ((X - cx).^2) ./ (rx^2) + ((Y - cy).^2) ./ (ry^2);
    opticDiskMask = normDist <= 1;
else
    opticDiskMask = false(origH, origW);
end

% Save image
ToolBox = getGlobalToolBox;
if ToolBox.getParams.saveFigures
    alpha = 0.2;
    color = [1, 0, 0];

    % Get ellipse
    t = linspace(0, 2*pi, 200);
    ex = bestBox(1) + (bestBox(3)/2) * cos(t);
    ey = bestBox(2) + (bestBox(4)/2) * sin(t);
    polyCoords = reshape([ex; ey], 1, []); % Interleave for insertShape
    outlineMask = insertShape(zeros(1024, 1024), 'Polygon', polyCoords, ...
        'Color', 'white', 'LineWidth', 4, 'Opacity', 1);
    isEdge = outlineMask(:,:,1) > 0.5;

    visImg = rgbM0; % This is 0-1 Single
    
    rCh = visImg(:,:,1);
    gCh = visImg(:,:,2);
    bCh = visImg(:,:,3);
    
    % Blend only where the edge exists
    rCh(isEdge) = rCh(isEdge) * (1 - alpha) + color(1) * alpha;
    gCh(isEdge) = gCh(isEdge) * (1 - alpha) + color(2) * alpha;
    bCh(isEdge) = bCh(isEdge) * (1 - alpha) + color(3) * alpha;
    
    visImg = cat(3, rCh, gCh, bCh);
    imwrite(visImg, fullfile(ToolBox.path_png, 'optic_disk.png'));

    maskImg = uint8(opticDiskMask) * 255;
    imwrite(maskImg, fullfile(ToolBox.path_png, 'optic_disk_mask.png'));
end

end