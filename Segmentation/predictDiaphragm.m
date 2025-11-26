function [diaphragm_mask, cx, cy, r] = predictDiaphragm(model, M0)
% predictDiaphragm - Predicts the diaphragm of a M0 image
%   Inputs:
%       model     - The pre-loaded YOLO model
%       M0        - Input image (2D matrix)
%   Output:
%       diaphragm_mask - Logical matrix (mask)
%       cx             - Diaphragm center x
%       cy             - Diaphragm center y
%       r              - Diaphragm radius

% Preprocess
inputSize = [1024, 1024]; 
[origH, origW] = size(M0);
resizedM0 = rescale(imresize(M0, inputSize));
rgbM0 = cat(3, resizedM0, resizedM0, resizedM0);
dlInput = dlarray(rgbM0, 'SSCB');

% Predict
[out1, out2] = predict(model, dlInput);

% Extract data
% out1: 1 x 37 x 21504 (Batch x Channels x Anchors)
% out2: 256 x 256 x 32 x 1    (H x W x MaskChannels x Batch)
dets = extractdata(gather(out1));
protos = extractdata(gather(out2));

% Generate Binary Mask
% This function handles all the YOLO decoding, sigmoid, and cropping
binaryMask = get_yolo_mask(dets, protos, inputSize);

% RANSAC circle
[contours, ~] = bwboundaries(binaryMask, 'noholes');

if ~isempty(contours)
    % Find largest contour
    [~, maxIdx] = max(cellfun(@length, contours));
    largestContour = contours{maxIdx}; 
    
    % Convert to x,y and run RANSAC
    pts = [largestContour(:,2), largestContour(:,1)];
    [cx, cy, r] = ransac_circle(pts, 500, 2.0);
    
end

% Output
if ~isempty(cx)
    [xx, yy] = meshgrid(1:inputSize(2), 1:inputSize(1));
    finalMask = ((xx - cx).^2 + (yy - cy).^2) <= r^2;
else
    finalMask = false(inputSize);
end

diaphragm_mask = imresize(finalMask, [origH, origW], 'nearest');

% Image save
ToolBox = getGlobalToolBox;
if ToolBox.getParams.saveFigures
    save_visualisation(M0, diaphragm_mask);
end

end

%% --- Helper Functions ---

function binaryMask = get_yolo_mask(dets, protos, targetSize)
% GET_YOLO_MASK Decodes YOLOv11 output into a binary mask
%   Inputs:
%       dets   - Detection head (1 x 37 x 21504)
%       protos - Prototype head (256 x 256 x 32 x 1)
%       targetSize - The size to resize the final mask to [1024 1024]

    % 1. Decode Detections
    dets = squeeze(dets); 
    
    % Find best anchor
    scores = dets(5, :); 
    [maxScore, idx] = max(scores);
    
    if maxScore < 0.25
        binaryMask = false(targetSize);
        return;
    end
    
    maskCoeffs = dets(6:37, idx); % 32x1
    box = dets(1:4, idx);         % [cx, cy, w, h]
    
    % 2. Decode Prototypes & Assemble
    protos = squeeze(protos); 
    % protos is now [256, 256, 32]
    
    % --- THE FIX: PERMUTE SPATIAL DIMS ---
    % Python is Row-Major, MATLAB is Column-Major. 
    % We swap the first two dimensions (Width and Height) to read memory correctly.
    protos = permute(protos, [2, 1, 3]); 
    
    [pw, ph, pc] = size(protos); % Note: dim 1 is now Width, dim 2 is Height
    
    % Flatten: Now traversing correctly
    protosFlat = reshape(protos, pw*ph, pc);
    
    % Matrix Multiplication
    maskRaw = protosFlat * maskCoeffs;
    
    % Reshape back
    maskImg = reshape(maskRaw, pw, ph);
    
    % Transpose back to normal (Height x Width) for MATLAB Image orientation
    maskImg = maskImg'; 
    
    % Sigmoid Activation
    maskImg = 1 ./ (1 + exp(-maskImg));
    
    % 3. Resize and Crop
    % Resize mask to 1024x1024
    maskFull = imresize(maskImg, targetSize, 'bilinear');
    
    % Crop to Bounding Box
    bx = box(1); by = box(2); bw = box(3); bh = box(4);
    
    x1 = max(1, round(bx - bw/2));
    y1 = max(1, round(by - bh/2));
    x2 = min(targetSize(2), round(bx + bw/2));
    y2 = min(targetSize(1), round(by + bh/2));
    
    boxMask = false(targetSize);
    boxMask(y1:y2, x1:x2) = true;
    
    % Final Binary Mask
    binaryMask = (maskFull > 0.5) & boxMask;
end

function [best_cx, best_cy, best_r] = ransac_circle(points, n_iter, tol)
    % Standard RANSAC logic
    num_points = size(points, 1);
    best_inliers = -1;
    best_cx = []; best_cy = []; best_r = [];

    if num_points < 3, return; end

    for i = 1:n_iter
        idx = randperm(num_points, 3);
        sample = points(idx, :);
        [cx, cy, r] = circle_from_points(sample);
        
        if isempty(cx) || isnan(cx) || r > 2000, continue; end
        
        dists = sqrt((points(:,1) - cx).^2 + (points(:,2) - cy).^2);
        inliers = sum(abs(dists - r) < tol);
        
        if inliers > best_inliers
            best_inliers = inliers;
            best_cx = cx; best_cy = cy; best_r = r;
        end
    end
end

function [cx, cy, r] = circle_from_points(p)
    x = p(:,1); y = p(:,2);
    x1 = x(1); x2 = x(2); x3 = x(3);
    y1 = y(1); y2 = y(2); y3 = y(3);
    
    temp = x2^2 + y2^2;
    bc = (x1^2 + y1^2 - temp) / 2;
    cd = (temp - x3^2 - y3^2) / 2;
    det = (x1 - x2) * (y2 - y3) - (x2 - x3) * (y1 - y2);
    
    if abs(det) < 1e-12
        cx = []; cy = []; r = []; return;
    end
    
    cx = (bc * (y2 - y3) - cd * (y1 - y2)) / det;
    cy = ((x1 - x2) * cd - (x2 - x3) * bc) / det;
    r = sqrt((x1 - cx)^2 + (y1 - cy)^2);
end

function save_visualisation(M0, diaphragm)
    alpha = 0.2;
    
    % 1. Prepare the base grayscale image (0 to 1 range)
    grayImg = rescale(im2double(M0));
    
    % 2. Convert to RGB by duplicating the gray channel 3 times
    % Size becomes [Height x Width x 3]
    outputImg = cat(3, grayImg, grayImg, grayImg);
    
    % 3. Define the Target Color (Green = [R=0, G=1, B=0])
    targetColor = [1, 0, 0]; 
    
    % 4. Loop through R, G, and B channels to apply the blend
    for c = 1:3
        channel = outputImg(:,:,c);
        
        channel(diaphragm) = (channel(diaphragm) * (1 - alpha)) + (targetColor(c) * alpha);
        
        % Update the output image channel
        outputImg(:,:,c) = channel;
    end

    % 5. Save
    ToolBox = getGlobalToolBox;
    imwrite(outputImg, fullfile(ToolBox.path_png, 'eye_diaphragm.png'));

    maskImg = uint8(diaphragm) * 255;
    imwrite(maskImg, fullfile(ToolBox.path_png, 'eye_diaphragm_mask.png'));
end