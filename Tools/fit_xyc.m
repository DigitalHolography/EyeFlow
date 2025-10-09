function [PWV, Tx, Ty, S, m, pks, idx, rows, cols] = fit_xyc(Z, dx, dy, name, branch_index)
% --- Input ---
% Z =  your 2D map (matrix: ny x nx)

ToolBox = getGlobalToolBox;

outputDir = fullfile(ToolBox.path_png, 'flexion');

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

if nargin < 2 dx = 1.0; end % pixel spacing in x (change if physical units known)

if nargin < 3 dy = 1.0; end % pixel spacing in y

    % --- Preprocess ---

    Z0 = Z - mean(Z(:)); % remove DC offset

    Z0(1:5, :) = [];
    Z0(end - 4:end, :) = [];

    [ny, nx] = size(Z0);
    cy = round(ny / 2) + 1;
    cx = round(nx / 2) + 1;

    Z0 = Z0 .* ~diskMask(ny, nx, 0.001); % remove central point
    % apply 2D Hanning (Hann) window to reduce edge effects before FFT
    wy = hann(ny);
    wx = hann(nx);
    W = wy * wx';
    Z0 = Z0 .* W;

    % --- FFT ---
    Nmult = 8; % must use for fine measurement
    F = fftshift(fft2(Z0,Nmult*ny,Nmult*nx));
    S = abs(F);

    bandWidth = diskMask(Nmult*ny, Nmult*nx, 0, 0.5);

    S = S .* bandWidth;

    [pks, rows, cols] = findpeaks2(S, 0, 0.8);

    % --- Frequency axes ---
    fx = (-Nmult*nx / 2:Nmult*nx / 2 - 1) / (Nmult*nx * dx); % cycles per unit in x
    fy = (-Nmult*ny / 2:Nmult*ny / 2 - 1) / (Nmult*ny * dy); % cycles per unit in y

    figure('Visible','on'), imagesc(fx,fy,S);
    axis xy;
    hold on;
    scatter(fx(cols), fy(rows), 80, 'o', 'r', 'LineWidth', 1.5);

    % --- Find dominant peak ---
    [m, idx] = max(pks);
    iy = rows(idx);
    ix = cols(idx);

    % --- Convert to radians per unit ---
    kx = 2 * pi * fx(ix); % a coefficient
    ky = 2 * pi * fy(iy); % b coefficient

    % Optional: compute spatial periods
    Tx = - 2 * pi / (kx); %% todo fix
    Ty = 2 * pi / (ky);
    fprintf('Period along x = %f units\n', Tx);
    fprintf('Period along y = %f units\n', Ty);

    % --- Plot map with wavevector direction ---
    figure('Visible', 'on');
    imagesc(linspace(-dx*nx/2, nx*dx/2,nx),linspace(-dy*ny/2, ny*dy/2,ny),Z0,[-0.1,0.1]);
    xlabel("time delay (s)");
    ylabel("arc length lag (mm)")
    hold on;

    % center point of the image
    xc = nx / 2;
    yc = ny / 2;

    if m > 10
        PWV = Ty / Tx;
        plot([0 - Tx / dx *0.2 , 0 + Tx / dx *0.2], [0 - Ty / dy *0.2, 0 + Ty / dy *0.2], 'w-', 'LineWidth', 2);
    else
        PWV = NaN;
    end

    text(0 , 0 , sprintf('%.2f mm / s', PWV), ...
        'Color', 'w', 'FontSize', 10, 'FontWeight', 'bold');

    
    
    title('Flexural Pulse Wave velocity measure');
    hold off;

    % Save figure
    saveas(gcf, fullfile(outputDir, ...
        sprintf("%s_%s_%d_signals_asym_over_time.png", ToolBox.folder_name, name, branch_index)));

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

    if ~isempty(pks) && minHeightRatio > 0
        M = max(pks);
        keep = pks > minHeightRatio * M;
        pks = pks(keep);
        rows = rows(keep);
        cols = cols(keep);
    end

end
