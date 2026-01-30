function Bout = imrotatecustom(A, angle)
% Rotates image A by a given angle, keeping output the same size,
% and filling undefined areas with NaN.

% Original size
sz = size(A);

% Create mask of defined (non-NaN) values
mask = ~isnan(A);

% Use a unique fill value that won't appear in A
% Handle case where all values are NaN
if any(mask(:))
    fill_value = min(A(mask)) - 1; % ensure it's below any real value
else
    fill_value = NaN; % fallback if entire image is NaN
end

% Replace NaNs with fill_value for rotation
A_temp = A;
A_temp(~mask) = fill_value;

% Rotate both the image and mask
Bfull = imrotate(A_temp, angle, 'bilinear', 'loose');
mask_rot = imrotate(double(mask), angle, 'bilinear', 'loose'); % Convert to double for rotation

% Set undefined regions to NaN
% To handle the cases where the A is a logical array (so no NaNs)
try
    Bfull(mask_rot < 0.5) = NaN;
catch
    Bfull(mask_rot < 0.5) = 0;
end

% Compute crop to center result and match original size
sz_new = size(Bfull);
r_start = max(1, round((sz_new(1) - sz(1)) / 2 + 1));
c_start = max(1, round((sz_new(2) - sz(2)) / 2 + 1));

% Ensure indices are within bounds
r_end = min(sz_new(1), r_start + sz(1) - 1);
c_end = min(sz_new(2), c_start + sz(2) - 1);

% Adjust if crop window exceeds rotated image size
if r_end > sz_new(1)
    r_end = sz_new(1);
    r_start = r_end - sz(1) + 1;
end

if c_end > sz_new(2)
    c_end = sz_new(2);
    c_start = c_end - sz(2) + 1;
end

% Crop to original size
Bout = Bfull(r_start:r_end, c_start:c_end);
end
