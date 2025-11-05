function mask = opticDiskDetect(model, img)
% opticDiskDetect - Detects the optic disk in a retinal image using a segmentation network
% Inputs:
%   model: The loaded eye diaphragm segmentation network
%   img: 2D matrix of the retinal image
% Outputs:
%   mask: Binary mask of the detected optic disk

pyenv; % Ensure Python environment is initialized

% Convert the image to a Python-compatible format
img_py = py.numpy.array(img);
% Call the Python function to perform segmentation
mask_py = py.yolo_segment_image(model, img_py);
% Convert the output mask back to MATLAB format
mask = double(mask_py);

end
