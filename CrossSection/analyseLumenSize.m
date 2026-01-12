function [fit_systole_cycles, fit_diastole_cycles] = analyseLumenSize(v_profile, systole_frames, diastole_frames, frame_gap_threshold)
% Takes v_profile in and return two arrays of arteries size

if nargin < 3
    error('Three inputs are required: v_profile, systole_frames, and diastole_frames.');
end

if nargin < 4
    frame_gap_threshold = 1;
end

ToolBox = getGlobalToolBox;
params = ToolBox.getParams;

systole_frame_groups = subdivide_array(systole_frames, frame_gap_threshold);
diastole_frame_groups = subdivide_array(diastole_frames, frame_gap_threshold);

systole_results = process_cycles(v_profile, systole_frame_groups);
diastole_results = process_cycles(v_profile, diastole_frame_groups);

get_diameter = @(s) get_diameter_or_nan(s);

systole_diameters = cellfun(get_diameter, systole_results);
diastole_diameters = cellfun(get_diameter, diastole_results);

fit_systole_cycles = systole_diameters * params.px_size;
fit_diastole_cycles = diastole_diameters * params.px_size;

end

function output_cell_array = subdivide_array(input_array, threshold)
% This function subdivides an input array into subarrays based on element spread.
% It groups consecutive numbers together.

if ~isvector(input_array)
    error('Input must be a vector.');
end

if nargin < 2
    threshold = 1;
end

d = diff(input_array);
break_points = find(d > threshold);

output_cell_array = cell(1, length(break_points) + 1);

if isempty(break_points)
    output_cell_array{1} = input_array;
    return;
end

start_index = 1;

for i = 1:length(break_points)
    end_index = break_points(i);
    output_cell_array{i} = input_array(start_index:end_index);
    start_index = end_index + 1;
end

output_cell_array{end} = input_array(start_index:end);
end

function diameter = get_diameter_or_nan(s)

if ~isempty(s)
    diameter = s.diameter;
else
    diameter = NaN;
end

end

function output_cell = process_cycles(v_profile, group_frames)
num_cycles = length(group_frames);
output_cell = cell(1, num_cycles);

for i = 1:num_cycles
    current_frames = group_frames{i};

    v_cycle = v_profile(:, current_frames);
    avg_profile = mean(v_cycle, 2, 'omitnan');

    try
        output_cell{i} = fit_parabol_diam(avg_profile);

        if output_cell{i}.rsquared < 0.8
            % TODO: Remove Warning and put NaN
            warning('R-squared value (%.2f) is low for cycle %d. The fit may not be reliable.', output_cell{i}.rsquared, i);
        end

    catch ME
        warning(ME.identifier, 'Failed to fit cycle profile for cycle %d. Error: %s', i, ME.message);
        output_cell{i} = []; % Store empty if fit fails
    end

end

end

% TODO: maybe make
function fit_results = fit_parabol_diam(velocity_profile)
% Fits a parabolic profile and calculates the diameter from its roots.
%
% This function is a wrapper for customPoly2Fit and customPoly2Roots.
% It takes a 1D velocity profile, fits a quadratic polynomial (parabola)
% to it, finds the roots of the polynomial, and calculates the distance
% between them, which represents the diameter.
%
% INPUT:
%   velocity_profile - A 1D column or row vector of velocity data.
%
% OUTPUT:
%   fit_results      - A struct containing all results:
%       .diameter                   : The calculated diameter (abs(root1 - root2)).
%       .diameter_error             : The propagated error of the diameter.
%       .rsquared                   : The R-squared value of the parabolic fit.
%       .p1, .p2, .p3               : Coefficients of the fitted polynomial.
%       .p1_err, .p2_err, .p3_err   : Errors of the coefficients.
%       .root1, .root2              : The two roots of the polynomial.
%       .root1_err, .root2_err      : The errors of the roots.
%

if isempty(velocity_profile) || numel(velocity_profile) < 3
    warning('Input velocity profile is too short for a parabolic fit. Returning empty.');
    fit_results = struct('diameter', NaN, 'diameter_error', NaN, 'rsquared', NaN);
    return;
end

% Ensure column vector
y = velocity_profile(:);

central_range = find(y > 0.1 * max(y));
centt = mean(central_range);

r_range = (central_range - centt);

[p1, p2, p3, rsquared, p1_err, p2_err, p3_err] = customPoly2Fit(r_range', y(central_range));

% Check poor fit
% if rsquared < 0.8
%     warning('R-squared value (%.2f) is low. The fit may not be reliable.', rsquared);
% end

[r1, r2, r1_err, r2_err] = customPoly2Roots(p1, p2, p3, p1_err, p2_err, p3_err);

diameter = abs(r1 - r2);
diameter_error = sqrt(r1_err ^ 2 + r2_err ^ 2);

fit_results = struct();

fit_results.diameter = diameter;
fit_results.diameter_error = diameter_error;
fit_results.rsquared = rsquared;
fit_results.p1 = p1;
fit_results.p2 = p2;
fit_results.p3 = p3;
fit_results.p1_err = p1_err;
fit_results.p2_err = p2_err;
fit_results.p3_err = p3_err;
fit_results.root1 = r1;
fit_results.root2 = r2;
fit_results.root1_err = r1_err;
fit_results.root2_err = r2_err;
fit_results.center_offset = centt;
fit_results.central_range = central_range;
end

% +=====================================================================+ %
% |                           DEBUG FUNCTIONS                           | %
% +=====================================================================+ %

% Debug function to show parabole and its fit
function show_para_and_fit(value, fit_results)
x_data = 1:length(value);

p = [fit_results.p1, fit_results.p2, fit_results.p3];

% Create a smooth x-axis for the fitted curve
x_fit = linspace(fit_results.root1, fit_results.root2, 200);

x_fit_original = x_fit + fit_results.center_offset;

% Calculate the corresponding y-values for the fitted curve
y_fit = polyval(p, x_fit);

figure;
hold on;

plot(x_data, value, '.', 'Color', [0.7 0.7 0.7], 'DisplayName', 'Full Profile Data');
plot(fit_results.central_range, value(fit_results.central_range), 'b.', 'MarkerSize', 12, 'DisplayName', 'Data Used for Fit');
plot(x_fit_original, y_fit, 'r-', 'LineWidth', 2, 'DisplayName', 'Parabolic Fit');

plot([fit_results.root1, fit_results.root2] + fit_results.center_offset, [0, 0], 'kx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Roots');

title('Comparison of Average Profile and Parabolic Fit');
xlabel('Spatial Points (e.g., Pixel Index)');
ylabel('Average Velocity (or other metric)');
legend('show', 'Location', 'best'); % Add a legend
grid on; % Add a grid for easier reading
hold off; % Release the plot
end
