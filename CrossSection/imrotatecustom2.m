function Bout = imrotatecustom2(A, angle)

sz = size(A);

mask = ~isnan(A);

if any(mask(:))
    fill_value = min(A(mask)) - 1;
else
    fill_value = NaN;
end

A_temp = A;
A_temp(~mask) = fill_value;

Bfull = imrotate(A_temp, angle, 'bilinear', 'loose');
mask_rot = imrotate(double(mask), angle, 'bilinear', 'loose') > 0.5;

Bfull(~mask_rot) = NaN;

rows = any(mask_rot, 2);
cols = any(mask_rot, 1);

r1 = find(rows, 1, 'first');
r2 = find(rows, 1, 'last');
c1 = find(cols, 1, 'first');
c2 = find(cols, 1, 'last');

Bcrop = Bfull(r1:r2, c1:c2);

out = NaN(sz);

rs = size(Bcrop,1);
cs = size(Bcrop,2);

r_start = floor((sz(1) - rs)/2) + 1;
c_start = floor((sz(2) - cs)/2) + 1;

r_idx = r_start:(r_start+rs-1);
c_idx = c_start:(c_start+cs-1);

valid_r = r_idx >= 1 & r_idx <= sz(1);
valid_c = c_idx >= 1 & c_idx <= sz(2);

out(r_idx(valid_r), c_idx(valid_c)) = ...
    Bcrop(valid_r, valid_c);

Bout = out;

end
