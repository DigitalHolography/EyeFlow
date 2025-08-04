function [D, dD, A, dA, c1, c2, rsquare] = computeVesselCrossSection(subImg, figName, ToolBox, papillaDiameter, flagFigures)

arguments
    subImg
    figName
    ToolBox
    papillaDiameter
    flagFigures = true
end

% Parameters
params = ToolBox.getParams;
HydrodynamicDiameters = params.json.CrossSectionsAnalysis.HydrodynamicDiameters;

[numX, ~] = size(subImg);

if ~isnan(papillaDiameter) && ~isempty(papillaDiameter)
    px_size = 1.8 / papillaDiameter / (2 ^ params.json.Preprocess.InterpolationFactor);
else
    px_size = params.px_size;
end

if ~HydrodynamicDiameters
    D = mean(sum(~isnan(subImg), 2)); % in pixels
    dD = 0;
    A = pi * (D * px_size / 2) ^ 2;
    dA = 0;
    c1 = 1;
    c2 = numX;
    rsquare = 1;
    return
end

% Compute velocity profile
profile = mean(subImg, 1, 'omitnan');
L = length(profile);

% Find all points above 50% threshold
central_range = find(profile > 0.1 * max(profile));
centt = mean(central_range);

r_range = (central_range - centt) * px_size;

[p1, p2, p3, rsquare, p1_err, p2_err, p3_err] = customPoly2Fit(r_range', profile(central_range)');
[r1, r2, r1_err, r2_err] = customPoly2Roots(p1, p2, p3, p1_err, p2_err, p3_err);

if r1 > r2
    r1 = NaN;
    r2 = NaN;
end

% Calculate cross-section limits in pixel indices
c1 = max(ceil(centt + (r1 / px_size)), 1);
c2 = min(floor(centt + (r2 / px_size)), L);

% Determine cross-section width
D = abs(r1 - r2);
dD = sqrt(r1_err ^ 2 + r2_err ^ 2);

if (D > sqrt(2) * L * px_size) || (rsquare < 0.6)
    D = NaN;
    dD = NaN;
end

% Compute cross-sectional area
A = pi * (D / 2) ^ 2;
dA = pi * (D / 2) .* dD;

if flagFigures
    poiseuilleProfileFigure(subImg, px_size, p1, p2, p3, r1, r2, figName, ToolBox);
end

end
