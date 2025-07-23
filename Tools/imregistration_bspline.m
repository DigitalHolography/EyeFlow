function [registered, ctrlPtsDisp] = imregistration_bspline(fixed, moving, gridSpacing, maxIter)
% IMREGISTRATION_BSPLINE Register moving to fixed image using B-spline FFD
%
% Inputs:
%   fixed       - fixed grayscale image (MxN)
%   moving      - moving grayscale image (MxN)
%   gridSpacing - control point spacing in pixels (e.g. 20)
%   maxIter     - max number of optimization iterations (e.g. 50)
%
% Outputs:
%   registered  - warped moving image aligned to fixed
%   ctrlPtsDisp - displacement vectors at B-spline control points [PxQx2]

% Convert images to double
fixed = im2double(fixed);
moving = im2double(moving);
[M, N] = size(fixed);

% Create B-spline control point grid
cx = ceil(M / gridSpacing) + 3; % +3 for cubic spline boundary
cy = ceil(N / gridSpacing) + 3;
ctrlPtsX = linspace(-gridSpacing, M + gridSpacing, cx);
ctrlPtsY = linspace(-gridSpacing, N + gridSpacing, cy);

% Initialize control point displacements to zero
ctrlPtsDisp = zeros(cx, cy, 2);

% Create B-spline basis knots for x and y
knotsX = augknt([ctrlPtsX(4:end - 3)], 3);
knotsY = augknt([ctrlPtsY(4:end - 3)], 3);

% Create meshgrid for image pixels
[X, Y] = meshgrid(1:N, 1:M);

% Optimization options
lr = 1; % learning rate for gradient descent

for iter = 1:maxIter
    % 1. Compute displacement field by evaluating B-spline at pixel locations
    dispX = bsplineEval(ctrlPtsDisp(:, :, 1), knotsX, knotsY, X, Y);
    dispY = bsplineEval(ctrlPtsDisp(:, :, 2), knotsX, knotsY, X, Y);

    % 2. Warp moving image with current displacement (bilinear interp)
    warpedMoving = interp2(moving, X + dispX, Y + dispY, 'linear', 0);

    % 3. Compute difference and gradient wrt displacement (simple SSD)
    diff = warpedMoving - fixed;

    % Compute image gradients of warped moving image
    [Ix, Iy] = gradient(warpedMoving);

    % 4. Compute gradient of cost wrt displacement field
    gradX = diff .* Ix;
    gradY = diff .* Iy;

    % 5. Backpropagate gradient to control points (chain rule)
    gradCtrlPtsX = bsplineBackprop(ctrlPtsDisp(:, :, 1), knotsX, knotsY, X, Y, gradX);
    gradCtrlPtsY = bsplineBackprop(ctrlPtsDisp(:, :, 2), knotsX, knotsY, X, Y, gradY);

    % 6. Update control points displacements (gradient descent)
    ctrlPtsDisp(:, :, 1) = ctrlPtsDisp(:, :, 1) - lr * gradCtrlPtsX;
    ctrlPtsDisp(:, :, 2) = ctrlPtsDisp(:, :, 2) - lr * gradCtrlPtsY;

    % Display iteration info
    if mod(iter, 10) == 0 || iter == 1
        loss = sum(diff(:) .^ 2);
        fprintf('Iter %d, SSD Loss: %.6f\n', iter, loss);
    end

end

% Final warp
dispX = bsplineEval(ctrlPtsDisp(:, :, 1), knotsX, knotsY, X, Y);
dispY = bsplineEval(ctrlPtsDisp(:, :, 2), knotsX, knotsY, X, Y);
registered = interp2(moving, X + dispX, Y + dispY, 'linear', 0);

end

%% Helper: Evaluate B-spline displacement field at pixel locations
function F = bsplineEval(ctrlPts, knotsX, knotsY, X, Y)
% ctrlPts: control points (cx x cy)
% knotsX, knotsY: knots vectors for B-spline basis
% X,Y: pixel coordinates

cx = size(ctrlPts, 1);
cy = size(ctrlPts, 2);

% Evaluate B-spline basis in x and y directions at pixel coordinates
Bx = bsplineBasis(knotsX, cx, X);
By = bsplineBasis(knotsY, cy, Y);

% Combine to get displacement field: F = Bx * ctrlPts * By'
F = zeros(size(X));

for i = 1:cx

    for j = 1:cy
        F = F + ctrlPts(i, j) * (Bx(:, i) .* By(:, j));
    end

end

F = reshape(F, size(X));
end

%% Helper: Compute B-spline basis functions for all points
function B = bsplineBasis(knots, nCtrlPts, coord)
% knots: knots vector
% nCtrlPts: number of control points
% coord: evaluation points (NxM)

% Using Cox-De Boor recursion to compute basis
% For efficiency, vectorize for 1D coord vector

coord = coord(:);
N = length(coord);
degree = 3; % cubic B-spline

B = zeros(N, nCtrlPts);

% Find knot span for each coord
for i = 1:nCtrlPts
    B(:, i) = bsplineBasisFunc(degree, i, knots, coord);
end

end

%% Helper: B-spline basis function (recursive Cox-de Boor)
function val = bsplineBasisFunc(p, i, knots, x)
% p: degree
% i: basis function index
% knots: knot vector
% x: points to evaluate (vector)

if p == 0
    val = double(knots(i) <= x & x < knots(i + 1));
else
    denom1 = knots(i + p) - knots(i);
    denom2 = knots(i + p + 1) - knots(i + 1);

    term1 = zeros(size(x));
    term2 = zeros(size(x));

    if denom1 > 0
        term1 = ((x - knots(i)) / denom1) .* bsplineBasisFunc(p - 1, i, knots, x);
    end

    if denom2 > 0
        term2 = ((knots(i + p + 1) - x) / denom2) .* bsplineBasisFunc(p - 1, i + 1, knots, x);
    end

    val = term1 + term2;
end

% Fix last knot interval
if i == length(knots) - p - 1
    val(x == knots(end)) = 1;
end

end

%% Helper: Backpropagate gradient to control points
function gradCtrlPts = bsplineBackprop(ctrlPts, knotsX, knotsY, X, Y, gradF)
% ctrlPts: current control point displacements
% gradF: gradient of loss wrt displacement field at pixels

cx = size(ctrlPts, 1);
cy = size(ctrlPts, 2);

Bx = bsplineBasis(knotsX, cx, X);
By = bsplineBasis(knotsY, cy, Y);

gradCtrlPts = zeros(cx, cy);

% Sum gradient contributions from all pixels weighted by basis functions
for i = 1:cx

    for j = 1:cy
        gradCtrlPts(i, j) = sum(sum(gradF .* (Bx(:, i) .* By(:, j))));
    end

end

end
