function VideoResizing(obj, params)
% Function to resize the video data cube based on specified parameters
% Inputs:
%   obj - PreprocessorClass object containing video data
%   params - Struct containing resizing parameters

out_height = params.json.Preprocess.Resize.FrameHeight;
out_width = params.json.Preprocess.Resize.FrameWidth;
out_numFrames = params.json.Preprocess.Resize.VideoLength;

[numX, numY, numFrames] = size(obj.M0_ff);

isHeightResized = out_height > 0;
isWidthResized = out_width > 0;
isLengthResized = out_numFrames > 0;

if isHeightResized && ~isWidthResized
    out_width = numX;

elseif ~isHeightResized && isWidthResized
    out_height = numY;

elseif ~isHeightResized && ~isWidthResized
    % If neither Height nor Width is inputed then the function makes the video
    % spatially isomorphic by resizing the smallest dimension to the largest
    out_width = max(numX, numY);
    out_height = max(numX, numY);

    if (numX == numY) && ~isLengthResized
        % If the video is spatially isomorphic and no Heigth/Width/Length
        % is inputed then the function has no more job to do
        return
    end

end

if ~isLengthResized
    out_numFrames = numFrames;
end

fprintf('Resizing data cube : %dx%dx%d to %dx%dx%d\n', ...
    numX, numY, numFrames, out_height, out_width, out_numFrames);

[Xq, Yq, Zq] = meshgrid(linspace(1, numY, out_width), linspace(1, numX, out_height), linspace(1, numFrames, out_numFrames));
% tmp_ref = zeros(numX, numY, numFrames);

tmp_calc = obj.M0;
tmp = interp3(tmp_calc, Xq, Yq, Zq);
obj.M0 = single(tmp);

tmp_calc = obj.M0_ff;
tmp = interp3(tmp_calc, Xq, Yq, Zq);
obj.M0_ff = single(tmp);

tmp_calc = obj.f_AVG;
tmp = interp3(tmp_calc, Xq, Yq, Zq);
obj.f_AVG = single(tmp);

tmp_calc = obj.f_RMS;
tmp = interp3(tmp_calc, Xq, Yq, Zq);
obj.f_RMS = single(tmp);
clear tmp_calc tmp Xq Yq Zq
end
