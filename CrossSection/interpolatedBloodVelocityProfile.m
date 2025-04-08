function interpolatedBloodVelocityProfile(v_cell, dv_cell, sysIdx, diasIdx, numSections, name)
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

Color_err = [0.7 0.7 0.7];

if strcmp(name, 'Artery')
    Color_sys = [1 0 0];
    Color_dias = [1/2 0 0];
else
    Color_sys = [0 0 1];
    Color_dias = [0 0 1/2];
end

% Get sizes
numCircles = params.json.CrossSectionsAnalysis.NumberOfCircles;
numFrames = size(v_cell{1}, 2);
numInterp = params.json.CrossSectionsFigures.InterpolationPoints;

w2w = linspace(-1, 1, numInterp);
x_label = 'wall-to-wall distance (a.u.)';
y_label = 'Velocity (mm/s)';

% Preallocate interpolated profiles
v_interp = cell(1, numCircles);
dv_interp = cell(1, numCircles);

% Interpolate profiles for all circles, sections, and frames
parfor cIdx = 1:numCircles
    v_interp{cIdx} = zeros(numSections(cIdx), numFrames, numInterp); % Preallocate for mean velocity
    dv_interp{cIdx} = zeros(numSections(cIdx), numFrames, numInterp); % Preallocate for uncertainties

    for sectionIdx = 1:numSections(cIdx)

        for frameIdx = 1:numFrames
            v = v_cell{cIdx}{sectionIdx, frameIdx}; % Mean velocity profile
            dv = dv_cell{cIdx}{sectionIdx, frameIdx}; % Velocity uncertainty profile

            % Handle edge case with negative velocities
            if any(v < 0)
                [~, locs] = findpeaks(-v); % Find peaks in negative velocities

                if isempty(locs)
                    indx = find(v > 0); % Use positive velocities if no peaks
                elseif length(locs) > 1
                    indx = locs(1):locs(end); % Use range between first and last peak
                else

                    if locs(1) > length(v) / 2
                        indx = 1:locs(1);
                    else
                        indx = locs(1):length(v);
                    end

                end

            else
                indx = find(v > 0); % Use positive velocities
            end

            % Interpolate profiles
            L = length(indx);
            v_interp{cIdx}(sectionIdx, frameIdx, :) = interp1(1:L, v(indx), linspace(1, L, numInterp)); % Interpolate mean velocity
            dv_interp{cIdx}(sectionIdx, frameIdx, :) = interp1(1:L, dv(indx), linspace(1, L, numInterp)); % Interpolate uncertainty
        end

    end

end

% Preallocate variables for systolic and diastolic profiles
v_sys = zeros(numInterp, 1);
dv_sys = zeros(numInterp, 1);
v_dias = zeros(numInterp, 1);
dv_dias = zeros(numInterp, 1);
v_video = zeros(numFrames, numInterp);
dv_video = zeros(numFrames, numInterp);

for cIdx = 1:numCircles

    for sectionIdx = 1:numSections(cIdx)
        % Average profiles for systolic and diastolic frames
        v_sys = v_sys + squeeze(sum(v_interp{cIdx}(sectionIdx, sysIdx, :), 2) / numSections(cIdx));
        dv_sys = dv_sys + squeeze(sum(dv_interp{cIdx}(sectionIdx, sysIdx, :) .^ 2, 2) / numSections(cIdx));
        v_dias = v_dias + squeeze(sum(v_interp{cIdx}(sectionIdx, diasIdx, :), 2) / numSections(cIdx));
        dv_dias = dv_dias + squeeze(sum(dv_interp{cIdx}(sectionIdx, diasIdx, :) .^ 2, 2) / numSections(cIdx));
        v_video = v_video + squeeze(v_interp{cIdx}(sectionIdx, :, :) / numSections(cIdx));
        dv_video = dv_video + squeeze(dv_interp{cIdx}(sectionIdx, :, :) .^ 2 / numSections(cIdx));
    end

end

numSys = length(sysIdx);
numDias = length(diasIdx);
v_sys = (v_sys / (numSys * numCircles))';
dv_sys = (sqrt(dv_sys) / (numSys * numCircles))';
v_dias = (v_dias / (numDias * numCircles))';
dv_dias = (sqrt(dv_dias) / (numDias * numCircles))';
v_video = (v_video / numCircles)';
dv_video = (sqrt(dv_video) / numCircles)';

% Create curves for plotting
curve1_sys = v_sys + dv_sys;
curve2_sys = v_sys - dv_sys;
inBetween_sys = [curve1_sys, flip(curve2_sys)];
ft2_sys = [w2w, flip(w2w)];

curve1_dias = v_dias + dv_dias;
curve2_dias = v_dias - dv_dias;
inBetween_dias = [curve1_dias, flip(curve2_dias)];
ft2_dias = [w2w, flip(w2w)];

% Plot systolic and diastolic profiles
figure("Visible", "off");
pbaspect([1.618 1 1]);
ax = gca;
hold on;

fill(ax, ft2_sys, inBetween_sys, Color_err, 'EdgeColor', 'none');
plot(ax, w2w, curve1_sys, 'Color', Color_err, 'LineWidth', 2);
plot(ax, w2w, curve2_sys, 'Color', Color_err, 'LineWidth', 2);
plot(ax, w2w, v_sys, '--', 'Color', Color_sys, 'LineWidth', 2);

fill(ax, ft2_dias, inBetween_dias, Color_err, 'EdgeColor', 'none');
plot(ax, w2w, curve1_dias, 'Color', Color_err, 'LineWidth', 2);
plot(ax, w2w, curve2_dias, 'Color', Color_err, 'LineWidth', 2);
plot(ax, w2w, v_dias, '--', 'Color', Color_dias, 'LineWidth', 2);

% Poiseuille fit
warning("off");
[~, centt_sys] = max(v_sys);
r_range_sys = (1:numInterp) - centt_sys;
f_sys = fit(r_range_sys', v_sys', 'poly2');
poiseuille_fit_sys = f_sys.p1 * r_range_sys .^ 2 + f_sys.p2 * r_range_sys + f_sys.p3;

[~, centt_dias] = max(v_dias);
r_range_dias = (1:numInterp) - centt_dias;
f_dias = fit(r_range_dias', v_dias', 'poly2');
poiseuille_fit_dias = f_dias.p1 * r_range_dias .^ 2 + f_dias.p2 * r_range_dias + f_dias.p3;

plot(ax, w2w, poiseuille_fit_sys, 'Color', Color_sys, 'LineWidth', 2);
plot(ax, w2w, poiseuille_fit_dias, 'Color', Color_dias, 'LineWidth', 2);
warning("on");

% Adjust axes and labels
axis padded;
axP = axis;
axis tight;
axT = axis;
axis([axT(1), axT(2), axP(3), 1.07 * axP(4)]);
hold off;

box on;
set(gca, 'Linewidth', 2);
xlabel('Profile (a.u.)', 'FontSize', 14);
ylabel('Velocity (mm/s)', 'FontSize', 14);
pbaspect([1.618 1 1]);

% Export figure
exportgraphics(gca, fullfile(ToolBox.path_png, 'crossSectionsAnalysis', sprintf("%s_diasys_%s.png", ToolBox.main_foldername, name)));

% Interpolated blood velocity profile and video export
if exportVideos
    fig = figure("Visible", "off");
    ax = axes(fig);
    hold(ax, "on");

    fillPlot = fill(ax, NaN, NaN, Color_err, 'EdgeColor', 'none');
    curve1Plot = plot(ax, NaN, NaN, "Color", Color_err, 'LineWidth', 2);
    curve2Plot = plot(ax, NaN, NaN, "Color", Color_err, 'LineWidth', 2);
    meanPlot = plot(ax, NaN, NaN, '--', "Color", Color_sys, 'LineWidth', 2);
    poiseuillePlot = plot(ax, NaN, NaN, '-', "Color", Color_sys, 'LineWidth', 2);
    xlabel(x_label);
    ylabel(y_label);

    pbaspect([1.618 1 1]);
    axis tight;
    ax.YLim = [-5; 30];
    box on;
    set(gca, 'Linewidth', 2);

    video = zeros(420, 560, 3, numFrames, 'single'); % Preallocate video array

    parfor frameIdx = 1:numFrames
        % Compute mean and RMS profiles
        curve1_video = (v_video(:, frameIdx) + dv_video(:, frameIdx))';
        curve2_video = (v_video(:, frameIdx) - dv_video(:, frameIdx))';
        inBetween_video = [curve1_video, flip(curve2_video)];
        ft2_video = [w2w, flip(w2w)];

        % Update plot elements
        set(fillPlot, 'XData', ft2_video, 'YData', inBetween_video);
        set(curve1Plot, 'XData', w2w, 'YData', curve1_video);
        set(curve2Plot, 'XData', w2w, 'YData', curve2_video);
        set(meanPlot, 'XData', w2w, 'YData', v_video(:, frameIdx));

        % Poiseuille fit
        warning("off");
        [~, centt] = max(v_video(:, frameIdx));
        r_range = (1:numInterp) - centt;
        f = fit(r_range', v_video(:, frameIdx), 'poly2');
        poiseuille_fit = f.p1 * r_range .^ 2 + f.p2 * r_range + f.p3;
        set(poiseuillePlot, 'XData', w2w, 'YData', poiseuille_fit);
        warning("on");
        % Capture frame
        video(:, :, :, frameIdx) = rescale(frame2im(getframe(fig)));
    end

    % Write video to disk
    writeGifOnDisc(video, sprintf("interp_profile_%s", name), "ToolBox", ToolBox);

end

end
