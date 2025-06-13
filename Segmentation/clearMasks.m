function [mask_dilated, mask_opened, mask_closed, mask_widened] = clearMasks(input_mask, params)
% clearMasks - Processes a binary mask to remove small objects, close gaps, and dilate.
%
% Inputs:
%   input_mask  - A binary mask (logical array) to be processed.
%   params      - A structure containing parameters for processing:
%       - min_area: Minimum area for objects to be retained.
%       - imclose_radius: Radius for morphological closing.
%       - min_width: Minimum width for vessels.
%       - imdilate_size: Size for final dilation.
% Outputs:
%   mask_dilated - The processed binary mask after morphological operations.
%   mask_opened   - The mask after area opening.
%   mask_closed   - The mask after morphological closing.
%   mask_widened  - The mask after ensuring minimum width.

arguments
    input_mask logical % Binary mask to be processed
    params.min_area = 50 % Minimum area for objects to be retained
    params.imclose_radius = 2 % Radius for morphological closing
    params.min_width = 1 % Minimum width for vessels
    params.imdilate_size = 2 % Size for final dilation
end

% Extract parameters from the input structure
min_area = params.min_area;
imclose_radius = params.imclose_radius;
min_width = params.min_width;
imdilate_size = params.imdilate_size;

% Step 1: Remove small objects using area opening
mask_opened = bwareaopen(input_mask, min_area, 4);

% Step 2: Close small gaps using morphological closing
imcloseSE = strel('disk', imclose_radius);
mask_closed = imclose(mask_opened, imcloseSE);

% Step 3: Ensure minimum mask width using skeletonization and dilation
minWidthSE = strel('disk', min_width);
skel = bwskel(mask_closed);
mask_widened = mask_closed | imdilate(skel, minWidthSE);

% Step 4: Final dilation
dilationSE = strel('disk', imdilate_size);
mask_dilated = imdilate(mask_widened, dilationSE);

end
