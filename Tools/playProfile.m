function playProfile(v_profile,x)
% INPUT:
%   v_profile - Array (numInterp x numFrames) containing the measured
% velocity profile data across x (freq or time).
% Get global toolbox settings
ToolBox = getGlobalToolBox;
params = ToolBox.getParams;
numFrames = size(v_profile, 2);
numInterp = params.json.CrossSectionsFigures.InterpolationPoints;
w2w = linspace(-1, 1, numInterp);

% Create figure for static plot
figure("Visible", "on");
hold('on');
ylim([min(real(v_profile),[],"all")-3,max(real(v_profile),[],"all")+3]);

% Plot profile data
imag_line = plot(w2w, imag(v_profile(:,1)), '-', 'Color', 'r', 'LineWidth', 2);
real_line = plot(w2w, real(v_profile(:,1)), '-', 'Color', 'b', 'LineWidth', 2);

for i = 1:numFrames
    set(imag_line, 'YData', imag(v_profile(:,i)));
    set(real_line, 'YData', real(v_profile(:,i)));
    pause(0.10)
end
end