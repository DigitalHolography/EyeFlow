function VideoInterpolating(obj)
% This function interpolates the video frames of M0_ff, f_AVG, and f_RMS
% using the specified interpolation factor from the parameters JSON file.
% It uses parfor for parallel processing of each frame to speed up the interpolation process.

[numX, numY, numFrames] = size(obj.M0_ff);
params = Parameters_json(obj.directory, obj.param_name);

kInterp = params.json.Preprocess.InterpolationFactor;

out_numX = (numX - 1) * (2 ^ kInterp - 1) + numX;
out_numY = (numY - 1) * (2 ^ kInterp - 1) + numY;

fprintf("        - Interpolating data cube : %dx%dx%d > %dx%dx%d\n", ...
    numX, numY, numFrames, out_numX, out_numY, numFrames);

% Reference M0
tmp_M0_ff = zeros(out_numX, out_numY, numFrames);
tmp_Calc_M0_ff = obj.M0_ff;

parfor frameIdx = 1:numFrames
    tmp_M0_ff(:, :, frameIdx) = interp2(tmp_Calc_M0_ff(:, :, frameIdx), kInterp);
end

obj.M0_ff = tmp_M0_ff;
clear tmp_M0_ff tmp_Calc_M0_ff

% M1M0
tmpM1M0 = zeros(out_numX, out_numY, numFrames);
tmpCalcM1M0 = obj.f_AVG;

parfor frameIdx = 1:numFrames
    tmpM1M0(:, :, frameIdx) = interp2(tmpCalcM1M0(:, :, frameIdx), kInterp);
end

obj.f_AVG = tmpM1M0;
clear tmpM1M0 tmpCalcM1M0

% M2M0
tmpM2M0 = zeros(out_numX, out_numY, numFrames);
tmpCalcM2M0 = obj.f_RMS;

parfor frameIdx = 1:numFrames
    tmpM2M0(:, :, frameIdx) = interp2(tmpCalcM2M0(:, :, frameIdx), kInterp);
end

obj.f_RMS = single(tmpM2M0);
clear tmpM2M0 tmpCalcM2M0

end
