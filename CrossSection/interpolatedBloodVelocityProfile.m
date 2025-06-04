function interpolatedBloodVelocityProfile(v_cell, dv_cell, sysIdx, diasIdx, name)
% interpolatedBloodVelocityProfile - Combines and optimizes the computation and plotting of velocity profiles.
% Inputs:
%   v_cell: Cell array containing mean velocity profiles for each circle, section, and frame.
%   dv_cell: Cell array containing velocity profile uncertainties.
%   sysIdx: Indices of systolic frames.
%   diasIdx: Indices of diastolic frames.
%   numSections: Number of sections for each circle.
%   name: Name identifier for saving files.
%   rad: Radius values for each circle.

% Get global toolbox settings
ToolBox = getGlobalToolBox;
params = ToolBox.getParams;
exportVideos = params.exportVideos;

% Set visualization parameters
Color_err = [0.7 0.7 0.7];

if strcmp(name, 'Artery')
    Color_sys = [1 0 0];
    Color_dias = [0.5 0 0];
else
    Color_sys = [0 0 1];
    Color_dias = [0 0 0.5];
end

% Get dimensions
[numCircles, numBranches] = size(v_cell);
numFrames = size(v_cell{1}, 2);
numInterp = params.json.CrossSectionsFigures.InterpolationPoints;
w2w = linspace(-1, 1, numInterp);

% Preallocate interpolated arrays
v_interp = cell(numCircles, numBranches);
vse_interp = cell(numCircles, numBranches);

% Process each circle and branch
for cIdx = 1:numCircles

    for bIdx = 1:numBranches

        if ~isempty(v_cell{cIdx, bIdx})
            current_v = v_cell{cIdx, bIdx};
            current_dv = dv_cell{cIdx, bIdx};

            % Preallocate for this branch
            temp_v = zeros(numFrames, numInterp);
            temp_dv = zeros(numFrames, numInterp);

            for frameIdx = 1:numFrames
                v = current_v{frameIdx};
                dv = current_dv{frameIdx};

                % Find valid indices (positive velocities)
                indx = [];

                if any(v < 0)
                    [~, locs] = findpeaks(-v);

                    if isempty(locs)
                        indx = find(v > 0);
                    elseif numel(locs) > 1
                        indx = locs(1):locs(end);
                    end

                else
                    indx = find(v > 0);
                end

                if isempty(indx)
                    temp_v(frameIdx, :) = zeros(1, numInterp);
                    temp_dv(frameIdx, :) = zeros(1, numInterp);
                elseif numel(indx) == 1
                    temp_v(frameIdx, :) = zeros(1, numInterp);
                    temp_dv(frameIdx, :) = zeros(1, numInterp);
                else
                    % Vectorized interpolation
                    L = numel(indx);
                    xq = linspace(1, L, numInterp);
                    temp_v(frameIdx, :) = interp1(1:L, v(indx), xq, 'linear');
                    temp_dv(frameIdx, :) = interp1(1:L, dv(indx), xq, 'linear');
                end

            end

            v_interp{cIdx, bIdx} = temp_v;
            vse_interp{cIdx, bIdx} = temp_dv;
        end

    end

end

% Initialize accumulation variables
v_sys_acc = zeros(1, numInterp);
dv_sys_acc = zeros(1, numInterp);
v_dias_acc = zeros(1, numInterp);
dv_dias_acc = zeros(1, numInterp);
v_video_acc = zeros(numFrames, numInterp);
dv_video_acc = zeros(numFrames, numInterp);

% Accumulate profiles
validBranches = 0;

for cIdx = 1:numCircles

    for bIdx = 1:numBranches

        if ~isempty(v_interp{cIdx, bIdx})
            validBranches = validBranches + 1;

            % Systolic accumulation
            v_sys_acc = v_sys_acc + sum(v_interp{cIdx, bIdx}(sysIdx, :), 1);
            dv_sys_acc = dv_sys_acc + sum(vse_interp{cIdx, bIdx}(sysIdx, :) .^ 2, 1);

            % Diastolic accumulation
            v_dias_acc = v_dias_acc + sum(v_interp{cIdx, bIdx}(diasIdx, :), 1);
            dv_dias_acc = dv_dias_acc + sum(vse_interp{cIdx, bIdx}(diasIdx, :) .^ 2, 1);

            % Video accumulation
            v_video_acc = v_video_acc + v_interp{cIdx, bIdx};
            dv_video_acc = dv_video_acc + vse_interp{cIdx, bIdx} .^ 2;
        end

    end

end

% Check if any valid branches were found

if validBranches == 0
    warning('No valid branches found for the specified systolic and diastolic indices.');
    return;
end

% Compute average profiles
numSys = numel(sysIdx);
numDias = numel(diasIdx);
v_sys = v_sys_acc / (numSys * validBranches);
dv_sys = sqrt(dv_sys_acc) / (numSys * validBranches);
v_dias = v_dias_acc / (numDias * validBranches);
dv_dias = sqrt(dv_dias_acc) / (numDias * validBranches);
v_video = v_video_acc / validBranches;
dv_video = sqrt(dv_video_acc) / validBranches;

% Create confidence bounds
createBounds = @(v, dv) struct( ...
    'upper', v + dv, ...
    'lower', v - dv, ...
    'x', [w2w, fliplr(w2w)], ...
    'y', [v + dv, fliplr(v - dv)]);

bounds_sys = createBounds(v_sys, dv_sys);
bounds_dias = createBounds(v_dias, dv_dias);

% Create figure for static plot
figure("Visible", "off");
hold('on');

% Plot systolic data
fill(bounds_sys.x, bounds_sys.y, Color_err, 'EdgeColor', 'none');
plot(w2w, bounds_sys.upper, 'Color', Color_err, 'LineWidth', 2);
plot(w2w, bounds_sys.lower, 'Color', Color_err, 'LineWidth', 2);
plot(w2w, v_sys, '--', 'Color', Color_sys, 'LineWidth', 2);

% Plot diastolic data
fill(bounds_dias.x, bounds_dias.y, Color_err, 'EdgeColor', 'none');
plot(w2w, bounds_dias.upper, 'Color', Color_err, 'LineWidth', 2);
plot(w2w, bounds_dias.lower, 'Color', Color_err, 'LineWidth', 2);
plot(w2w, v_dias, '--', 'Color', Color_dias, 'LineWidth', 2);

% Add Poiseuille fits
warning('off', 'curvefit:fit:noStartPoint');
fitPoiseuille = @(v) fit((1:numInterp)' - (numInterp / 2), v', 'poly2');
f_sys = fitPoiseuille(v_sys);
f_dias = fitPoiseuille(v_dias);

r_range = (1:numInterp) - (numInterp / 2); % Center around midpoint
plot(w2w, f_sys.p1 * r_range .^ 2 + f_sys.p2 * r_range + f_sys.p3, ...
    'Color', Color_sys, 'LineWidth', 2);
plot(w2w, f_dias.p1 * r_range .^ 2 + f_dias.p2 * r_range + f_dias.p3, ...
    'Color', Color_dias, 'LineWidth', 2);
warning('on', 'curvefit:fit:noStartPoint');

% Finalize static plot
axis('tight');
ylim([min([bounds_sys.lower, bounds_dias.lower]), ...
               1.07 * max([bounds_sys.upper, bounds_dias.upper])]);
xlabel('wall-to-wall distance (a.u.)', 'FontSize', 14);
ylabel('Velocity (mm/s)', 'FontSize', 14);
pbaspect([1.618 1 1]);
box('on');
set(gca, 'LineWidth', 2);

% Export static figure
outputDir = fullfile(ToolBox.path_png, 'crossSectionsAnalysis');

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

exportgraphics(gca, fullfile(outputDir, ...
    sprintf("%s_diasys_%s.png", ToolBox.folder_name, name)), 'Resolution', 300);

% Video export if requested
if exportVideos
    figVideo = figure('Visible', 'off');
    hold('on');

    % Create plot elements
    fillPlot = fill(NaN, NaN, Color_err, 'EdgeColor', 'none');
    upperPlot = plot(NaN, NaN, 'Color', Color_err, 'LineWidth', 2);
    lowerPlot = plot(NaN, NaN, 'Color', Color_err, 'LineWidth', 2);
    meanPlot = plot(NaN, NaN, '--', 'Color', Color_sys, 'LineWidth', 2);
    fitPlot = plot(NaN, NaN, '-', 'Color', Color_sys, 'LineWidth', 2);

    % Configure axes
    axis('tight');
    ylim(get(gca, 'YLim'));
    xlabel('wall-to-wall distance (a.u.)');
    ylabel('Velocity (mm/s)');
    pbaspect([1.618 1 1]);
    box('on');
    set(gca, 'LineWidth', 2);

    % Preallocate video
    % Create bounds for this frame
    bounds_frame = createBounds(v_video(1, :), dv_video(1, :));

    % Update plots
    set(fillPlot, 'XData', bounds_frame.x, 'YData', bounds_frame.y);
    set(upperPlot, 'XData', w2w, 'YData', bounds_frame.upper);
    set(lowerPlot, 'XData', w2w, 'YData', bounds_frame.lower);
    set(meanPlot, 'XData', w2w, 'YData', v_video(1, :));

    ylim([min([bounds_sys.lower, bounds_dias.lower]), ...
                   1.07 * max([bounds_sys.upper, bounds_dias.upper])]);

    [numX_fig, numY_fig, ~] = size(frame2im(getframe(figVideo)));
    video = zeros(numX_fig, numY_fig, 3, numFrames);

    % Process each frame
    for frameIdx = 1:numFrames
        % Create bounds for this frame
        bounds_frame = createBounds(v_video(frameIdx, :), dv_video(frameIdx, :));

        % Update plots
        set(fillPlot, 'XData', bounds_frame.x, 'YData', bounds_frame.y);
        set(upperPlot, 'XData', w2w, 'YData', bounds_frame.upper);
        set(lowerPlot, 'XData', w2w, 'YData', bounds_frame.lower);
        set(meanPlot, 'XData', w2w, 'YData', v_video(frameIdx, :));

        % Update Poiseuille fit
        f_frame = fitPoiseuille(v_video(frameIdx, :));
        set(fitPlot, 'XData', w2w, ...
            'YData', f_frame.p1 * r_range .^ 2 + f_frame.p2 * r_range + f_frame.p3);

        % Capture frame
        video(:, :, :, frameIdx) = frame2im(getframe(figVideo));
    end

    % Write video
    writeGifOnDisc(mat2gray(video), sprintf("wall2wall_profile_%s", name), "ToolBox", ToolBox);
end

end
