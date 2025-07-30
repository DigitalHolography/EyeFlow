function Bout = imrotatecustom(A, angle)
    % Rotates image A by a given angle, keeping output the same size,
    % and filling undefined areas with NaN.

    % Original size
    sz = size(A);

    % Create mask of defined (non-NaN) values
    mask = ~isnan(A);

    % Use a unique fill value that won't appear in A
    fill_value = min(A(~isnan(A))) - 1;  % ensure it's below any real value

    % Replace NaNs with fill_value for rotation
    A_temp = A;
    A_temp(~mask) = fill_value;

    % Rotate both the image and mask
    Bfull = imrotate(A_temp, angle, 'bilinear', 'loose');
    mask_rot = imrotate(mask, angle, 'bilinear', 'loose');

    % Set undefined regions to NaN
    Bfull(mask_rot < 0.5) = NaN;

    % Compute crop to center result and match original size
    sz_new = size(Bfull);
    r_start = round((sz_new(1) - sz(1)) / 2) + 1;
    c_start = round((sz_new(2) - sz(2)) / 2) + 1;

    % Ensure indices are within bounds
    r_end = r_start + sz(1) - 1;
    c_end = c_start + sz(2) - 1;

    if r_start < 1 || c_start < 1 || r_end > sz_new(1) || c_end > sz_new(2)
        error('Rotated image cannot be cropped to original size. Try reducing angle.');
    end

    % Crop to original size
    Bout = Bfull(r_start:r_end, c_start:c_end);
end
