function [PWV, Tx, Ty, S, m, pks, idx, rows, cols] = fit_xyc(Z, dx, dy, name, branch_index)
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

    %% --- Get global toolbox and prepare output directory ---
    ToolBox = getGlobalToolBox();
    outputDir = fullfile(ToolBox.path_png, 'flexion');
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

    %% --- Preprocess ---
    Z0 = Z - mean(Z(:), 'omitnan'); % remove DC offset and NaN safety
    if size(Z0,1) <= 10 || size(Z0,2) <= 10
        error('Input map is too small.');
    end

    % Remove edges (optional cleaning)
    Z0(1:5, :) = [];
    Z0(end-4:end, :) = [];

    [ny, nx] = size(Z0);

    % Mask out the central region if needed
    Z0 = Z0 .* ~diskMask(ny, nx, 0.001);
    wy = hann(ny);
    wx = hann(nx);
    W = wy * wx';
    Z0 = Z0 .* W;

    %% --- FFT ---
    Nmult = 8; % zero-padding multiplier for higher frequency resolution
    F = fftshift(fft2(Z0, Nmult*ny, Nmult*nx));
    S = abs(F);

    % Keep only central frequencies
    bandWidth = diskMask(Nmult*ny, Nmult*nx, 0, 0.5);
    S = S .* bandWidth;
    S0 = S;
    S0(1:ny*Nmult > ny*Nmult/2,:) = 0;
    S0(:,1:nx*Nmult < nx*Nmult/2) = 0;

    %% --- Peak detection ---
    [pks, rows, cols] = findpeaks2(S0, 0, 0.25);
    if isempty(pks)
        %warning('No peaks found in spectrum. Returning NaN results.');
        PWV = NaN; Tx = NaN; Ty = NaN; m = NaN;
        idx = []; rows = []; cols = [];
        return;
    end

    %% --- Frequency axes ---
    fx = (-Nmult*nx/2 : Nmult*nx/2-1) / (Nmult*nx*dx); % cycles per unit x
    fy = (-Nmult*ny/2 : Nmult*ny/2-1) / (Nmult*ny*dy); % cycles per unit y

    %% --- Plot frequency spectrum ---
    fig1 = figure('Visible', 'on');
    imagesc(fx, fy, S);
    axis xy; colormap turbo; colorbar;
    hold on;
    scatter(fx(cols), fy(rows), 80, 'ro', 'LineWidth', 1.5);
    xlabel('Spatial frequency (1/x)');
    ylabel('Spatial frequency (1/y)');
    title('2D FFT Spectrum and Peak Detection');

    %% --- Dominant peak ---
    [m, idx] = max(pks);
    iy = rows(idx);
    ix = cols(idx);

    kx = 2 * pi * fx(ix); % radians per unit
    ky = 2 * pi * fy(iy);

    % Avoid division by zero
    if abs(kx) < eps || abs(ky) < eps
        warning('Invalid frequency detected (zero or near zero).');
        Tx = NaN; Ty = NaN; PWV = NaN;
        return;
    end

    %% --- Compute periods and PWV ---
    Tx = 2 * pi / abs(kx);
    Ty = 2 * pi / abs(ky);

    % fprintf('Period along x = %.3f units\n', Tx);
    % fprintf('Period along y = %.3f units\n', Ty);

    PWV = Ty / Tx;

    if PWV < 1 || PWV > 4 % mm/s
        PWV = NaN;
    end

    %% --- Plot original map with wave direction ---
    fig2 = figure('Visible', 'on');
    xVals = linspace(-dx*nx/2, dx*nx/2, nx);
    yVals = linspace(-dy*ny/2, dy*ny/2, ny);
    imagesc(xVals, yVals, Z0, [-0.1, 0.1]);
    axis xy; colormap parula; colorbar;
    xlabel('time delay (s)');
    ylabel('arc length lag (mm)');
    hold on;

    % Draw wavevector direction
    quiver(0, 0, Tx*0.2, Ty*0.2, 0, 'w', 'LineWidth', 2, 'MaxHeadSize', 2);
    text(0, 0, sprintf('%.2f mm/s', PWV), 'Color', 'w', 'FontSize', 10, 'FontWeight', 'bold');
    title('Flexural Pulse Wave Velocity');

    %% --- Save figures ---
    saveas(fig1, fullfile(outputDir, sprintf('%s_%s_%d_fft_spectrum.png', ToolBox.folder_name, name, branch_index)));
    saveas(fig2, fullfile(outputDir, sprintf('%s_%s_%d_wave_map.png', ToolBox.folder_name, name, branch_index)));

    %close(fig1); close(fig2);
end


%% --- Helper: find peaks in 2D ---
function [pks, rows, cols] = findpeaks2(Z, smoothSigma, minHeightRatio)
    if nargin < 2, smoothSigma = 0; end
    if nargin < 3, minHeightRatio = 0; end

    if smoothSigma > 0
        Z = imgaussfilt(Z, smoothSigma);
    end

    BW = imregionalmax(Z);
    [rows, cols] = find(BW);
    pks = Z(BW);

    if isempty(pks)
        return;
    end

    if minHeightRatio > 0
        M = max(pks);
        keep = pks > minHeightRatio * M;
        pks = pks(keep);
        rows = rows(keep);
        cols = cols(keep);
    end
end
