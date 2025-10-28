function [alphaWom, pseudoViscosity, fitParams] = WomersleyNumberEstimation(v_profile, cardiac_frequency, name, idx, circleIdx, branchIdx)
% WomersleyNumberEstimation estimates the dimensionless Womersley number (alphaWom)
% by fitting the velocity profile to a Womersley flow model.
%
% This improved version fits a 5-parameter model to find the optimal:
% 1. alpha: Womersley number
% 2. amplitude: Complex scaling factor (magnitude and phase)
% 3. center: Center position of the vessel
% 4. width: Effective radius of the vessel
%
% INPUT:
%   v_profile         - Array (crossSectionLength x numFrames) of velocity data.
%   cardiac_frequency - Cardiac frequency in Hz.
%   name, idx, ...    - Identifiers for saving the output plot.
%
% OUTPUT:
%   alphaWom          - Estimated Womersley number.
%   pseudoViscosity   - Derived dynamic viscosity in Pa·s.
%   fitParams         - A struct containing all fitted parameters.


ToolBox = getGlobalToolBox;
params = ToolBox.getParams;
NUM_INTERP_POINTS = params.json.exportCrossSectionResults.InterpolationPoints;
FFT_PADDING_FACTOR = 16;

alphaWom = NaN;
pseudoViscosity = NaN;
fitParams = struct('alpha', NaN, 'amplitude', NaN, 'center', NaN, 'width', NaN);

v_profile_avg = mean(v_profile, 2);
valid_idxs = v_profile_avg ~= 0;
v_profile = v_profile(valid_idxs, :);
crossSectionLength = size(v_profile, 1);

if crossSectionLength > 1
    v_profile = interp1(linspace(1, crossSectionLength, crossSectionLength), v_profile, linspace(1, crossSectionLength, NUM_INTERP_POINTS));
else
    warning('Not enough valid points in the velocity profile. Skipping fit.');
    return;
end

numFrames = size(v_profile, 2);
N_fft = numFrames * FFT_PADDING_FACTOR;
v_profile_ft = fftshift(fft(v_profile, N_fft, 2), 2);

f = fft_freq_vector(ToolBox.fs * 1000 / ToolBox.stride, N_fft);

[~, cardiac_idx] = min(abs(f - cardiac_frequency));
% Average over a small frequency band around the cardiac frequency for stability
freq_resolution = f(2) - f(1);
margin_hz = freq_resolution * 1.5; % Capture main lobe of the peak
cardiac_idxs = find(abs(f - cardiac_frequency) <= margin_hz);

if isempty(cardiac_idxs)
    warning('Cardiac frequency not found in FFT spectrum. Using closest peak.');
    cardiac_idxs = cardiac_idx;
end

% FIT PART


v_meas = mean(v_profile_ft(:, cardiac_idxs), 2);


x_coords = linspace(-1, 1, NUM_INTERP_POINTS);


uWom_base = @(alpha, r) (1 - (besselj(0, 1i^(3/2) * alpha * r) ./ besselj(0, 1i^(3/2) * alpha)));

costFun = @(p) [real(generate_womersley_model(p, x_coords, uWom_base) - v_meas.'); ...
                imag(generate_womersley_model(p, x_coords, uWom_base) - v_meas.')];


alpha_init = 4;
amp_init_complex = mean(v_meas(abs(v_meas)>0));
center_init = 0;
width_init = 0.8;

p_init = [alpha_init, real(amp_init_complex), imag(amp_init_complex), center_init, width_init];

lb = [0.1,  -Inf, -Inf, -0.8, 0.1]; % alpha > 0, width > 0
ub = [20,   Inf,  Inf,  0.8, 1.5];  % Reasonable physiological/physical bounds

options = optimoptions('lsqnonlin', 'Display', 'off', 'Algorithm', 'trust-region-reflective');
try
    [p_fit, ~] = lsqnonlin(costFun, p_init, lb, ub, options);

    alphaWom = p_fit(1);
    fitParams.alpha = p_fit(1);
    fitParams.amplitude = p_fit(2) + 1i * p_fit(3);
    fitParams.center = p_fit(4);
    fitParams.width = p_fit(5);

    omega = 2 * pi * cardiac_frequency;
    pseudoViscosity = NaN;

catch ME
    warning('Womersley fit failed for %s (idx %d): %s', name, idx, ME.message);
    % Return NaNs if fit fails
    return;
end

% --- Plotting and Saving Results ---

uWom_fit = generate_womersley_model(p_fit, x_coords, uWom_base);

hFig = figure("Visible", "off");
hold on;
title(sprintf('Womersley Fit for %s (idx %d)', name, idx), 'Interpreter', 'none');
plot(x_coords, real(v_meas), 'b-', 'LineWidth', 2, 'DisplayName', 'Measured Data (Real)');
plot(x_coords, imag(v_meas), 'r-', 'LineWidth', 2, 'DisplayName', 'Measured Data (Imag)');
plot(x_coords, real(uWom_fit), 'b--', 'LineWidth', 2, 'DisplayName', 'Model Fit (Real)');
plot(x_coords, imag(uWom_fit), 'r--', 'LineWidth', 2, 'DisplayName', 'Model Fit (Imag)');
hold off;

xlim([-1 1]);
xlabel('Normalized Cross-section', 'FontSize', 14);
ylabel('Complex Velocity (a.u.)', 'FontSize', 14);
legend('show', 'Location', 'best');
box on;
grid on;
axis tight;
set(gca, 'LineWidth', 1.5);

% Add an annotation box with the key results
fit_string = sprintf('α Womersley: %.2f\nCenter: %.2f\nWidth: %.2f', ...
                     fitParams.alpha, fitParams.center, fitParams.width);
annotation('textbox', [0.2 0.7 0.2 0.2], 'String', fit_string, ...
            'FitBoxToText', 'on', 'BackgroundColor', 'w', ...
            'EdgeColor', 'k', 'FontSize', 12);

% Export static figure
save_path = fullfile(ToolBox.path_png, 'Womersley');
if ~isfolder(save_path)
    mkdir(save_path);
end
save_filename = fullfile(save_path, sprintf("%s_WomersleyFit_%s_idx%d_c%d_b%d.png", ToolBox.folder_name, name, idx, circleIdx, branchIdx));

try
    exportgraphics(hFig, save_filename, 'Resolution', 300);
catch export_error
    warning('Could not save figure');
end

% Close the figure if it was created invisibly
if ~strcmpi(get(hFig, 'Visible'), 'on')
    close(hFig);
end

end

    
function model_profile = generate_womersley_model(p, x, uWom_base)
% GENERATE_WOMERSLEY_MODEL Creates a scaled and positioned Womersley velocity profile.
%
% This function takes a parameter vector 'p' and generates the corresponding
% theoretical Womersley velocity profile on a given spatial coordinate axis 'x'.
%
% INPUTS:
%   p          - A 1x5 vector of model parameters:
%                p(1): alpha (Womersley number)
%                p(2): amplitude_real (real part of the complex scaling factor)
%                p(3): amplitude_imag (imaginary part of the complex scaling factor)
%                p(4): center (center position of the vessel on the x-axis)
%                p(5): width (effective radius of the vessel)
%
%   x          - A 1xN vector representing the fixed spatial coordinates where
%                the profile should be evaluated.
%
%   uWom_base  - A function handle to the base Womersley equation.
%                It should accept (alpha, r) as inputs, where 'r' is the
%                normalized radius from -1 to 1.
%
% OUTPUT:
%   model_profile - A 1xN complex-valued vector representing the generated
%                   Womersley velocity profile.

% 1. Unpack parameters for clarity
alpha         = p(1);
amplitude     = p(2) + 1i * p(3); % Reconstruct complex amplitude
center        = p(4);
width         = p(5);

% 2. Map the fixed 'x' coordinates to the vessel's normalized radius 'r'.
% This transforms the spatial coordinates into the model's coordinate system
% where the vessel exists from r = -1 to r = +1.
r = (x - center) / width;

% 3. Calculate the theoretical, unitless Womersley profile using the base function.
profile = uWom_base(alpha, r);

% 4. Enforce the no-slip boundary condition. The model is only valid
% inside the vessel walls (where abs(r) <= 1). Outside the walls, the
% velocity must be zero.
profile(abs(r) > 1) = 0;

% 5. Scale the unit profile by the complex amplitude to match the data's
% magnitude and overall phase. Assign to the output variable.
model_profile = amplitude * profile;

end

% =========================================================================
% OLD IMPLEM
% =========================================================================


% function [alphaWom, pseudoViscosity] = WomersleyNumberEstimation(v_profile, cardiac_frequency, name, idx, circleIdx, branchIdx)
% % WomersleyNumberEstimation estimates the dimensionless Womersley number (alphaWom)
% % by fitting the input velocity profile (v_profile) to a Womersley flow profile.
% %
% % INPUT:
% %   v_profile - Array (crossSectionLength x numFrames) containing the measured
% % velocity profile data across time.
% %   (dv - std)
% %   cardiac_frequency - Cardiac Frequency in Hz
% %
% % OUTPUT:
% %   alphaWom  - Estimated Womersley number (dimensionless), characterizing
% %               the pulsatile flow regime.
% %   pseudoViscosity  - induced dynamic viscosity
% % Create figure for static plot
% % Get global ToolBox settings
% ToolBox = getGlobalToolBox;
% params = ToolBox.getParams;
% crossSectionLength = size(v_profile,1);
% numFrames = size(v_profile, 2);
% 
% t = ToolBox.Cache.t;
% N_padding = 16;
% 
% alphaWom = NaN;
% pseudoViscosity = NaN;
% 
% numInterp = params.json.exportCrossSectionResults.InterpolationPoints;
% 
% 
% % v_profile = interp1(linspace(1,crossSectionLength,crossSectionLength),v_profile,linspace(1,crossSectionLength,numInterp));
% 
% % Trimming the 0 on the edges
% v_profile_avg  = squeeze(mean(v_profile, 2));
% idxs = v_profile_avg~=0;
% v_profile_avg = v_profile_avg(idxs);
% v_profile = v_profile(idxs, :);
% 
% % Interpolation for more points
% sectionSize = length(v_profile_avg);
% v_profile = interp1(linspace(1, sectionSize, sectionSize), v_profile, linspace(1, sectionSize, numInterp));
% 
% % Center with a weighted mean
% v_profile_avg  = squeeze(mean(v_profile, 2));
% indxCenter = sum((1:numInterp)' .* v_profile_avg) / sum(v_profile_avg);
% shiftIndx = indxCenter - (numInterp+1)/2;
% v_profile = circshift(v_profile, round(shiftIndx), 1);
% 
% v_profile_ft = fftshift(fft(v_profile, numFrames * N_padding, 2), 2);
% 
% w2w = linspace(-1, 1, numInterp);
% 
% f = fft_freq_vector(ToolBox.fs * 1000 / ToolBox.stride, numFrames * N_padding);
% 
% [~, cardiac_idx] = min(abs(f - cardiac_frequency));
% 
% margin_ = round(0.01 * numFrames * N_padding / 2 ); % +- 10 % fs /2
% cardiac_idxs = cardiac_idx + (-margin_:margin_);
% cardiac_idxs(cardiac_idxs > numFrames * N_padding) = [];
% cardiac_idxs(cardiac_idxs < 1) = [];
% 
% v_meas = mean(v_profile_ft(:, cardiac_idxs), 2);
% % Force the two Womersley hypothesis - bad hypothesis
% % v_profile_hyp = setSymetry(v_profile);
% 
% % v_profile_hyp = setBoundariesZeros(v_profile_hyp);  - bad hypothesis
% 
% % Calculate Fourier transform and display
% 
% % figure(1655), imagesc(t, w2w,v_profile); xlabel('Time'), ylabel('Vessel Cross Section');
% 
% % figure(1786), imagesc(f, w2w,log10(abs(v_profile_ft))); xlabel('Freq'), ylabel('Vessel Cross Section ');
% 
% % figure(1917), imagesc(angle((v_profile_ft))); xlabel('Freq'), ylabel('Vessel Cross Section ');
% 
% v_norm = v_meas / mean(v_meas);
% 
% % v_norm = v_norm - 2i * imag(v_norm);
% 
% % Create figure for static plot
% hFig = figure("Visible", "off");
% hold('on');
% 
% % Plot profile data
% plot(w2w, imag(v_norm), '-', 'Color', 'r', 'LineWidth', 2); 
% plot(w2w, real(v_norm), '-', 'Color', 'b', 'LineWidth', 2);
% 
% % Fit Womersley analytical profile
% R = 1;
% r = linspace(-1, 1, numInterp) *R; % normalized radius
% 
% uWom = @(alpha, r) (1 - (besselj(0, 1i ^ (3/2) * alpha * r) ./ besselj(0, 1i ^ (3/2) * alpha)));
% 
% % regulation_window = @(r, R) max(1 - (r / R) .^ 2, 0); % increase the weight of the central values in the fitting
% 
% uWom_tofit = @(alpha) (uWom(alpha, r) / mean(uWom(alpha, r)));
% 
% % costFun = @(alpha) norm(regulation_window(r, 1) .* uWom_tofit(alpha) - v_norm'); % least square error minimization
% costFun = @(alpha) norm(uWom_tofit(alpha) - v_norm'); % least square error minimization
% 
% alpha_init = 3; % initial guess
% alphaWom = fminsearch(costFun, alpha_init);
% % alphaWom = lsqnonlin(costFun, alpha_init, [], []);
% 
% uWom_fit = uWom_tofit(alphaWom);
% plot(w2w, imag(uWom_fit), '--', 'Color', 'r', 'LineWidth', 2);
% plot(w2w, real(uWom_fit), '--', 'Color', 'b', 'LineWidth', 2);
% 
% 
% 
% % plot(w2w, regulation_window(r, R), '--', 'Color', 'k', 'LineWidth', 2);
% 
% % Finalize static plot
% xlim([-1 1]);
% 
% xlabel('lumen cross-section (a.u.)', 'FontSize', 14);
% ylabel('Velocity (mm/s)', 'FontSize', 14);
% 
% box on
% axis tight
% axis padded;
% set(gca, 'LineWidth', 2);
% set(gca, 'PlotBoxAspectRatio', [1.618, 1, 1])
% ax = gca;
% ax.LineStyleOrderIndex = 1; % Reset if needed
% ax.SortMethod = 'depth'; % Try changing sorting method
% ax.Layer = 'top'; % This may help in some cases
% 
% annotation('textbox', [0.5 0.6 0.2 0.1], ...
%             'String', sprintf('alpha Womersley : %.1f  ± %.1f', alphaWom, NaN), ...
%             'FitBoxToText', 'on', ...
%             'BackgroundColor', 'w', ...
%             'EdgeColor', 'none', ...
%             'FontSize', 12);
% 
% % Export static figure
% 
% ax = gca;
% 
% if ~isfolder(fullfile(ToolBox.path_png, 'Womersley'))
%     mkdir(fullfile(ToolBox.path_png, 'Womersley'))
% end
% 
% if isvalid(ax)
%     exportgraphics(gca, fullfile(ToolBox.path_png, "Womersley", sprintf("%s_WomersleyFit_%s_idx%d_c%d_b%d.png", ToolBox.folder_name, name, idx, circleIdx, branchIdx)), 'Resolution', 300);
% else
%     warning('Current axes are not valid. Skipping export.');
% end
% 
% % Close the figure if not needed
% if ~strcmpi(get(hFig, 'Visible'), 'on')
%     close(hFig);
% end
% 
% end
