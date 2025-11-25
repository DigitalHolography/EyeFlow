function opticDiskMask = predictOpticDisk(model, M0)
% predictOpticDisk - Predicts the optic disk of a M0 image
%   Inputs:
%       model     - The pre-loaded YOLO model
%       M0        - Input image (2D matrix)
%   Output:
%       opticDiskMask - Logical matrix (mask)

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

end