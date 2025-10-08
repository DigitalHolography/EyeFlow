function [Tx, Ty] = fit_xyc(Z,dx,dy)
% --- Input ---
% Z =  your 2D map (matrix: ny x nx)

ToolBox = getGlobalToolBox;



if nargin < 2 dx = 1.0; end % pixel spacing in x (change if physical units known)
if nargin < 3 dy = 1.0; end% pixel spacing in y

% --- Preprocess ---


Z0 = Z - mean(Z(:)); % remove DC offset

Z0(1:5,:)=[];
Z0(end-4:end,:)=[];

[ny, nx] = size(Z0);
cy = round(ny / 2) + 1;
cx = round(nx / 2) + 1;

Z0 = Z0 .* ~diskMask(ny,nx,0.001); % remove central point
% win = hann(ny) * hann(nx)';  % 2D window
% Z0 = Z0 .* win;

% --- FFT ---
F = fftshift(fft2(Z0));
S = abs(F);


% 

bandWidth = diskMask(ny, nx, 0,0.5);

S = S .* bandWidth ;



[pks, rows, cols] = findpeaks2(S, 1, 0.8);

% --- Find dominant peak ---
[~, idx] = max(pks);
iy = rows(idx);
ix = cols(idx);

disp(ix)
disp(iy)

% --- Frequency axes ---
fx = (-nx / 2:nx / 2 - 1) / (nx * dx); % cycles per unit in x
fy = (-ny / 2:ny / 2 - 1) / (ny * dy); % cycles per unit in y

% --- Convert to radians per unit ---
kx = 2 * pi * fx(ix); % a coefficient
ky = 2 * pi * fy(iy); % b coefficient

% Optional: compute spatial periods
Tx = 2 * pi / (kx);
Ty = 2 * pi / (ky);
fprintf('Period along x = %f units\n', Tx);
fprintf('Period along y = %f units\n', Ty);

% --- Plot map with wavevector direction ---
figure('Visible', 'on');
imagesc(Z0, [-0.1 0.1]); axis image; colormap jet;
hold on;

% center point of the image
xc = nx / 2;
yc = ny / 2;

% direction vector from (a,b)
scale = min(nx, ny) / 4; % length scaling for visibility
vx = scale * sign(kx); % just show direction
vy = scale * sign(ky);

plot([xc - Tx / dx, xc + Tx / dx], [yc - Ty / dy, yc + Ty / dy], 'w-', 'LineWidth', 2);

text(xc + vx, yc + vy, sprintf('%.2f %s over %s', Ty/Tx, 'yunit', 'xunit'), ...
    'Color', 'w', 'FontSize', 10, 'FontWeight', 'bold');

title('2D Map with Wavevector Direction');
axis off;
hold off;

end

function [pks, rows, cols] = findpeaks2(Z, smoothSigma, minHeightRatio)
    if nargin < 2, smoothSigma = 0; end
    if nargin < 3, minHeightRatio = 0; end
    if smoothSigma > 0
        Z = imgaussfilt(Z, smoothSigma);
    end
    
    BW = imregionalmax(Z);
    [rows, cols] = find(BW);
    pks = Z(BW);

    figure('Visible', 'on');
    imagesc(Z); axis image; colormap jet;
    
    if ~isempty(pks) && minHeightRatio > 0
        M = max(pks);
        s_idx = ~(pks > minHeightRatio * M);
        pks(s_idx) = []; rows(s_idx) = []; cols(s_idx) = [];
    end



end
