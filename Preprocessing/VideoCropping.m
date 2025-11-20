function [firstFrame, lastFrame] = VideoCropping(obj, firstFrame, lastFrame)
%Crop a video (matrix dim 3)

[~, ~, numFrames] = size(obj.M0);

if lastFrame == -1
    lastFrame = numFrames;
end

if firstFrame == -1
    firstFrame = 1;
end

if firstFrame >= 1 && firstFrame <= numFrames ...
        || lastFrame >= 1 && lastFrame <= numFrames ...
        && firstFrame < lastFrame

    obj.M0 = obj.M0(:, :, firstFrame:lastFrame);
    obj.M1 = obj.M1(:, :, firstFrame:lastFrame);
    obj.M2 = obj.M2(:, :, firstFrame:lastFrame);

    if ~isempty(obj.SH)
        obj.SH = obj.SH(:, :, :, firstFrame:lastFrame);
    end

    fprintf('Data cube frame: %d/%d to %d/%d\n', firstFrame, numFrames, lastFrame, numFrames);

elseif firstFrame < 1 || firstFrame > numFrames
    error('First frame (%d) is out of bounds. It should be between 1 and %d.', firstFrame, numFrames);
elseif lastFrame < 1 || lastFrame > numFrames
    error('Last frame (%d) is out of bounds. It should be between 1 and %d.', lastFrame, numFrames);
elseif firstFrame >= lastFrame
    error('First frame (%d) should be less than last frame (%d).', firstFrame, lastFrame);
end

end
