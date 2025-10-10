function PWV = pulseWaveVelocity(U, mask, branch_index, name)
% Computes the pulse wave velocity based on a cross correlation computation
% U is the field over which we compute the velocity and mask is the mask of
% the selected retinal artery

% U(x,y,t) usually M0
% center the [x,y] barycenter (the center of the CRA)
ToolBox = getGlobalToolBox;
params = ToolBox.getParams;
x_bary = ToolBox.Cache.xy_barycenter(1);
y_bary = ToolBox.Cache.xy_barycenter(2);
[numX, numY] = size(mask);
N_frame = size(U, 3);
PWV = NaN; % initialize output

outputDir = fullfile(ToolBox.path_png, 'flexion');

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

%% create a grid of points to select points along the skeletton of the artery mask
%
dxx = 2;
skel = bwskel(mask);
grid = ones([numY, numX]) < 0;
grid(1:dxx:end, :) = true;
grid(:, 1:dxx:end) = true;
interpoints = grid & skel; % get points interpolating with grid

numpoints = sum(interpoints, 'all'); % get the number of points

if numpoints < 1 % no intersection with grid (typically the mask is too small)
    return
end

%% register the positions of points stating by the one closest to the CRA then going from closest to closest

[interpoints_y, interpoints_x] = ind2sub(size(interpoints), find(interpoints)); % y first
k = dsearchn([interpoints_x, interpoints_y], [x_bary, y_bary]); % get the nearest point to the center

absx = zeros([1, numpoints]); % x and y position register
absy = zeros([1, numpoints]);
abs_dist = zeros([1, numpoints]); % vessel curvilign absis

absx(1) = interpoints_x(k); % nearest point to the center
absy(1) = interpoints_y(k);

figure('Visible', 'off');
imagesc (interpoints)
scatter(x_bary, y_bary, 80, 'o', 'r', 'LineWidth', 1.5);
scatter(absx(1), absy(1), 80, 'o', 'g', 'LineWidth', 1.5);

interpoints_x(k) = [];
interpoints_y(k) = [];
abs_dist(1) = 0;

for kb = 2:numpoints
    k = dsearchn([interpoints_x, interpoints_y], [absx(kb - 1), absy(kb - 1)]);
    absx(kb) = interpoints_x(k);
    absy(kb) = interpoints_y(k);
    interpoints_x(k) = []; % deleting point from list
    interpoints_y(k) = [];
end

for kb = 2:numpoints
    abs_dist(kb) = abs_dist(kb - 1) + sqrt((absx(kb) - absx(kb - 1)) ^ 2 + (absy(kb) - absy(kb - 1)) ^ 2) * params.px_size;
end

%% for the positions extract a signal

% with orhtogonal sections
halfwidth = 15;

% edges = bwboundaries(mask);
% edges = edges{1};

% edges_mask = false(size(mask));
% edges_mask(sub2ind(size(mask), edges)) = true;

L = zeros(size(mask));
numinterp = 50;
U_x = zeros(numpoints, size(U, 3), numinterp);
Ux_edge = zeros(numpoints, size(U, 3));

prev_line = [];

for i = 2:numpoints - 1
    % tangent/normal
    tx = absx(i + 1) - absx(i - 1);
    ty = absy(i + 1) - absy(i - 1);
    if tx == 0 && ty == 0, continue; end
    tangent = [tx, ty] / norm([tx, ty]);
    normal = [-tangent(2), tangent(1)];

    % endpoints of current line
    P3 = [absx(i) - halfwidth * normal(1), absy(i) - halfwidth * normal(2)];
    P4 = [absx(i) + halfwidth * normal(1), absy(i) + halfwidth * normal(2)];

    if ~isempty(prev_line)
        P1 = prev_line(1, :); % previous start
        P2 = prev_line(2, :); % previous end

        % parameter grid (controls resolution of strip filling)
        nu = 2 * halfwidth + 1;
        nv = round(sqrt(sum((P3 - P1) .^ 2))); % distance between strips
        [u, v] = meshgrid(linspace(0, 1, nu), linspace(0, 1, nv));

        % bilinear interpolation of quadrilateral
        X = (1 - v) .* ((1 - u) * P1(1) + u * P2(1)) + v .* ((1 - u) * P3(1) + u * P4(1));
        Y = (1 - v) .* ((1 - u) * P1(2) + u * P2(2)) + v .* ((1 - u) * P3(2) + u * P4(2));

        % round to pixel indices
        X = round(X); Y = round(Y);

        % keep inside image/mask
        inside = X >= 1 & X <= size(mask, 2) & Y >= 1 & Y <= size(mask, 1);
        X = X(inside); Y = Y(inside);
        idx = sub2ind(size(mask), Y, X);
        idx = idx(mask(idx));

        if ~isempty(idx)
            L(idx) = i;

            for j = 1:N_frame
                idx_t = sub2ind(size(U), Y(:), X(:), repelem(j, length(Y))');
                U_x(i, j, :) = interp1((1:length(U(idx_t))) / length(U(idx_t)), U(idx_t), (1:numinterp) / numinterp);
                line = U(sub2ind(size(U), X(:), Y(:), repelem(j, length(Y), 1)));
                Ux_edge(i, j) = mean(squeeze(line') .* linspace(-1, 1, length(line)));
            end

        end

    end

    prev_line = [P3; P4]; % save endpoints for next step
end

figure('Visible', 'on');
imagesc(L)
title('selected sections along the artery')

% Save figure
saveas(gcf, fullfile(outputDir, ...
    sprintf("%s_%s_%d_orthogonal_overlay.png", ToolBox.folder_name, name, branch_index)));
axis off;

Ux = Ux_edge;
Ux_n = (Ux - mean(Ux, 2)) ./ std(Ux, [], 2);

figure('Visible', 'on');

imagesc(linspace(0, N_frame * (ToolBox.stride / ToolBox.fs / 1000), N_frame), linspace(0, abs_dist(end), numpoints), real(Ux_n));
xlabel("time (s)");
ylabel("arc length (mm)")
% Save figure
saveas(gcf, fullfile(outputDir, ...
    sprintf("%s_%s_%d_signals_asym_over_time.png", ToolBox.folder_name, name, branch_index)));

Ux_n = Ux_n';
Ux_n(isnan(Ux_n)) = 0;
R = xcorr(Ux_n, 'unbiased');

[Nlags, ~] = size(R);
M = size(Ux_n, 2);
Ravg = zeros(Nlags, 2 * M + 1);

for i = -M:M
    idx = find_indices_compact(M, i);
    Ravg(:, i + M + 1) = mean(real(R(:, idx)), 2);
end

Ravg(isnan(Ravg)) = 0;
Ravg = Ravg';
figure('Visible', 'off'), imagesc(Ravg, [-0.1 0.1]);
axis off;

fit_xyc(Ravg, (ToolBox.stride / ToolBox.fs / 1000), (abs_dist(end) / numpoints), name, branch_index);

close all;
end
