function [PWV, dPWV, Tx, Ty, S, m, pks, idx, rows, cols] = fit_xyc(Z, dx, dy, name, branch_index)
%FIT_XYC  Estimate spatial frequency components and compute Pulse Wave Velocity (PWV)
%
% --- Input ---
% Z              : 2D map (matrix: ny x nx)
% dx, dy         : Pixel spacing in x and y (default = 1)
% name           : Identifier string for output file
% branch_index   : Index of the branch analyzed
%
% --- Output ---
% PWV            : Pulse Wave Velocity estimate
% Tx, Ty         : Period along x and y
% S              : 2D FFT magnitude spectrum
% m              : Maximum peak value
% pks, idx       : Peaks and index of the dominant one
% rows, cols     : Coordinates of detected peaks

%% --- Default inputs ---
if nargin < 2 || isempty(dx), dx = 1.0; end
if nargin < 3 || isempty(dy), dy = 1.0; end
if nargin < 4, name = 'unnamed'; end
if nargin < 5, branch_index = 0; end

%% --- Get global ToolBox and prepare output directory ---
ToolBox = getGlobalToolBox();
outputDir = fullfile(ToolBox.path_png, 'flexion');
params = ToolBox.getParams;
saveFigures = params.saveFigures;

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% --- Preprocess ---
Z0 = Z - mean(Z(:), 'omitnan'); % remove DC offset and NaN safety

if size(Z0, 1) <= 10 || size(Z0, 2) <= 10
    % warning('Input map is too small for FFT analysis. Returning NaN results.');
    PWV = NaN; Tx = NaN; Ty = NaN; m = NaN;
    idx = []; rows = []; cols = [];
    S = []; pks = [];
    return
end

% Remove edges (optional cleaning)
Z0(1:5, :) = [];
Z0(end - 4:end, :) = [];

[ny, nx] = size(Z0);

% Mask out the central region if needed
Z0 = Z0 .* ~diskMask(ny, nx, 0.001);
wy = hann(ny);
wx = hann(nx);
W = wy * wx';
Z0 = Z0 .* W;

% --- FFT ---
Nmult = 8; % zero-padding multiplier for higher frequency resolution
F = fftshift(fft2(Z0, Nmult * ny, Nmult * nx));
S = abs(F);

% Keep only central frequencies
bandWidth = diskMask(Nmult * ny, Nmult * nx, 0, 0.5);
S = S .* bandWidth;
S0 = S;
S0(1:ny * Nmult > ny * Nmult / 2, :) = 0;
S0(:, 1:nx * Nmult < nx * Nmult / 2) = 0;

% --- Frequency axes ---
fx = (-Nmult * nx / 2:Nmult * nx / 2 - 1) / (Nmult * nx * dx); % cycles per unit x
fy = (-Nmult * ny / 2:Nmult * ny / 2 - 1) / (Nmult * ny * dy); % cycles per unit y

% --- Peak detection ---
[pks, rows, cols, fx_err, fy_err] = findpeaks2(S0, fx, fy, 0.8);

if isempty(pks)
    % warning('No peaks found in spectrum. Returning NaN results.');
    PWV = NaN; Tx = NaN; Ty = NaN; m = NaN;
    idx = []; rows = []; cols = [];
    fx_err = NaN; fy_err = NaN;
    return;
end

% --- Plot frequency spectrum ---
fig1 = figure('Visible', 'off');
h_wave_map = imagesc(fx, fy, S);
xlabel('(Hz)');
ylabel('(mm-1)');
axis xy equal tight;

xlim([-10, 10]);
ylim([-10, 10]);
colormap turbo;
colorbar;
hold on;

% --- Dominant peak ---
[m, idx] = max(pks);
iy = rows(idx);
ix = cols(idx);

scatter(fx(round(ix)), fy(round(iy)), 80, 'ro', 'LineWidth', 1.5, 'DisplayName', 'Detected Peaks');

fx_peak = fx(round(ix));
fy_peak = fy(round(iy));

% --- Estimated uncertainty on frequencies ---
dfx = fx_err(idx);
dfy = fy_err(idx);

% Convert to radians per unit
kx = 2 * pi * fx_peak;
ky = 2 * pi * fy_peak;
dkx = 2 * pi * dfx;
dky = 2 * pi * dfy;

% Avoid division by zero
if abs(kx) < eps || abs(ky) < eps
    % warning('Invalid frequency detected (zero or near zero).');
    Tx = NaN; Ty = NaN; PWV = NaN; dPWV = NaN;
    return;
end

% --- Compute periods and PWV ---
Tx = 2 * pi / abs(kx);
Ty = 2 * pi / abs(ky);

dTx = abs(Tx * (dkx / kx));
dTy = abs(Ty * (dky / ky));

% fprintf('Period along x = %.3f units\n', Tx);
% fprintf('Period along y = %.3f units\n', Ty);

PWV = Ty / Tx;
dPWV = PWV * sqrt((dTx / Tx) ^ 2 + (dTy / Ty) ^ 2);

ToolBox.Output.Extra.add(sprintf("PulseWaveVelocity/PWV%s_%d/PWV", name, branch_index),PWV);
ToolBox.Output.Extra.add(sprintf("PulseWaveVelocity/PWV%s_%d/PWV_std", name, branch_index),dPWV);
ToolBox.Output.Extra.add(sprintf("PulseWaveVelocity/PWV%s_%d/PWV_unit", name, branch_index),"mm/s");
% s_wave_map = saveImagescToStruct(h_wave_map);
% ToolBox.Output.Extra.add(sprintf("PulseWaveVelocity/PWV%s_%d/PWV_intercorr", name, branch_index),s_wave_map);


if PWV < 1 || PWV > 4 % mm/s
    PWV = NaN; dPWV = NaN;
end

fprintf('PWV = %.3f ± %.3f (mm/s)\n', PWV, dPWV);

if saveFigures
    % --- Plot original map with wave direction ---
    fig2 = figure('Visible', 'off');
    xVals = linspace(-dx * nx / 2, dx * nx / 2, nx);
    yVals = linspace(-dy * ny / 2, dy * ny / 2, ny);
    imagesc(xVals, yVals, Z0, [-0.1, 0.1]);
    axis xy; colormap parula; colorbar;
    xlabel('time delay (s)');
    ylabel('arc length lag (mm)');
    hold on;

    % Draw wavevector direction
    quiver(0, 0, Tx * 0.2, Ty * 0.2, 0, 'w', 'LineWidth', 2, 'MaxHeadSize', 2);
    text(0, 0, sprintf('%.2f mm/s', PWV), 'Color', 'w', 'FontSize', 10, 'FontWeight', 'bold');
    title('Flexural Pulse Wave Velocity');

    % --- Save figures ---
    saveas(fig1, fullfile(outputDir, sprintf('%s_%s_%d_fft_spectrum.png', ToolBox.folder_name, name, branch_index)));
    saveas(fig2, fullfile(outputDir, sprintf('%s_%s_%d_wave_map.png', ToolBox.folder_name, name, branch_index)));

    close(fig1); close(fig2);

end

end

%% --- Helper: find peaks in 2D ---
function [pks, rows, cols, fx_err, fy_err] = findpeaks2(S, fx, fy, minHeightRatio)
%FINDPEAKS2 Detect 2D peaks and estimate subpixel positions and uncertainties
%
% Inputs:
%   S               - 2D magnitude spectrum
%   fx, fy          - frequency axis vectors
%   minHeightRatio  - keep peaks > minHeightRatio * max(S)
%
% Outputs:
%   pks              - peak values
%   rows, cols       - subpixel peak positions (indices)
%   fx_err, fy_err   - uncertainty in fx, fy

if nargin < 4, minHeightRatio = 0; end

BW = imregionalmax(S);
[rows, cols] = find(BW);
pks = S(BW);

if isempty(pks)
    fx_err = []; fy_err = [];
    return;
end

if minHeightRatio > 0
    M = max(pks);
    keep = pks > minHeightRatio * M;
    pks = pks(keep);
    rows = rows(keep);
    cols = cols(keep);
end

% Initialize uncertainty arrays
fx_err = zeros(size(rows));
fy_err = zeros(size(rows));

% Subpixel refinement using local quadratic fitting
for k = 1:numel(rows)
    i = rows(k); j = cols(k);

    % Skip border peaks
    if i <= 1 || i >= size(S, 1) || j <= 1 || j >= size(S, 2)
        fx_err(k) = NaN; fy_err(k) = NaN;
        continue;
    end

    % Extract 3x3 neighborhood
    local = double(S(i - 1:i + 1, j - 1:j + 1));

    % Fit 2D paraboloid z = a*x^2 + b*y^2 + c*x*y + d*x + e*y + f
    [X, Y] = meshgrid(-1:1, -1:1);
    A = [X(:) .^ 2, Y(:) .^ 2, X(:) .* Y(:), X(:), Y(:), ones(9, 1)];
    coeff = A \ local(:);

    a = coeff(1); b = coeff(2); c = coeff(3); d = coeff(4); e = coeff(5);
    denom = 4 * a * b - c ^ 2;

    if abs(denom) < eps
        dx = 0; dy = 0;
    else
        dx = (c * e - 2 * b * d) / denom;
        dy = (c * d - 2 * a * e) / denom;
    end

    % Clamp offset to avoid runaway fits
    dx = max(min(dx, 1), -1);
    dy = max(min(dy, 1), -1);

    % Subpixel position
    rows(k) = i + dy;
    cols(k) = j + dx;

    % Simple error estimate: distance from integer grid
    % (could be refined with residuals of quadratic fit)
    fx_step = mean(diff(fx));
    fy_step = mean(diff(fy));
    fx_err(k) = abs(fx_step * 0.5); % ~half grid uncertainty
    fy_err(k) = abs(fy_step * 0.5);
end

end
