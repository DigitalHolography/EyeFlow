function [Tx, Ty] = fit_xyc(Z)
% --- Input ---
% Z =  your 2D map (matrix: ny x nx)

ToolBox = getGlobalToolBox;

[ny, nx] = size(Z);

dx = 1.0; % pixel spacing in x (change if physical units known)
dy = 1.0; % pixel spacing in y

% --- Preprocess ---
Z0 = Z - mean(Z(:)); % remove DC offset

% --- FFT ---
F = fftshift(fft2(Z0));
S = abs(F);

% --- Remove DC neighborhood ---
cy = round(ny / 2) + 1;
cx = round(nx / 2) + 1;
S(cy - 2:cy + 2, cx - 2:cx + 2) = 0;

% --- Find dominant peak ---
[~, idx] = max(S(:));
[iy, ix] = ind2sub(size(S), idx);

% --- Frequency axes ---
fx = (-nx / 2:nx / 2 - 1) / (nx * dx); % cycles per unit in x
fy = (-ny / 2:ny / 2 - 1) / (ny * dy); % cycles per unit in y

% --- Convert to radians per unit ---
kx = 2 * pi * fx(ix); % a coefficient
ky = 2 * pi * fy(iy); % b coefficient

fprintf('a (kx) = %f rad/unit\n', kx);
fprintf('b (ky) = %f rad/unit\n', ky);

% Optional: compute spatial periods
Tx = 2 * pi / abs(kx);
Ty = 2 * pi / abs(ky);
fprintf('Period along x = %f units\n', Tx);
fprintf('Period along y = %f units\n', Ty);

% --- Plot map with wavevector direction ---
figure('Visible', 'off');
imagesc(Z, [-0.1 0.1]); axis image; colormap jet;
hold on;

% center point of the image
xc = nx / 2;
yc = ny / 2;

% direction vector from (a,b)
scale = min(nx, ny) / 4; % length scaling for visibility
vx = scale * sign(kx); % just show direction
vy = scale * sign(ky);

plot([xc - Tx, xc + Tx], [yc - Ty, yc + Ty], 'w-', 'LineWidth', 2);

text(xc + vx, yc + vy, sprintf('%.2f rad/%s over %s', ky, 'yunit', 'xunit'), ...
    'Color', 'w', 'FontSize', 10, 'FontWeight', 'bold');

title('2D Map with Wavevector Direction');
axis off;
colormap(ToolBox.Cache.cmapArtery);

hold off;

end
