function VideoRegistering(obj, firstFrame, lastFrame)
% VideoRegistering - Registers the video using intensity based registration

video = obj.M0_ff;
numX = size(video, 1);
numY = size(video, 2);

disc_ratio = 0.35; % parametrize this coef if needed

disc = diskMask(numX, numY, disc_ratio);
video_reg = video .* disc - disc .* sum(video .* disc, [1, 2]) / nnz(disc); % minus the mean in the disc of each frame
video_reg = reshape(video_reg, size(video, 1), size(video, 2), 1, size(video, 3)); % insert a dimension to match reegistration functions

video_reg = video_reg ./ (max(abs(video_reg), [], [1, 2])); % rescaling each frame but keeps mean at zero

image_ref = mean(video_reg(:, :, firstFrame:lastFrame), 3); % ref image is from 10 to 20 for example
[~, shifts, scores] = register_video_from_reference(video_reg, image_ref);

f1 = figure("Visible", 'off'); plot(scores / mean(scores), 'k'); ylim([0 max(0.1, 1.2 * max(scores) / mean(scores))]);
title('Registration Correlation Score (1 is good) (u.a.)');
saveas(f1, fullfile(obj.directory, 'eyeflow', sprintf("%s_%s", obj.filenames, 'RegistrationCorrelationScore.png')));
f2 = figure("Visible", 'off');
plot(shifts(1, :)); hold on; plot(shifts(2, :));
imwrite(f2, fullfile(obj.directory, 'eyeflow', sprintf("%s_%s", obj.filenames, 'RegistrationShiftsXY.png')));
close all

obj.M0_ff = register_video_from_shifts(video, shifts);

obj.M0 = register_video_from_shifts(obj.M0, shifts);
obj.M1 = register_video_from_shifts(obj.M1, shifts);
obj.M2 = register_video_from_shifts(obj.M2, shifts);

end
