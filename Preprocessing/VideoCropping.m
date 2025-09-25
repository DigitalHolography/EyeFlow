function VideoCropping(obj)
%Crop a video (matrix dim 3)
params = Parameters_json(obj.directory, obj.param_name);
firstFrame = params.json.Preprocess.Crop.StartFrame;
lastFrame = params.json.Preprocess.Crop.EndFrame;

if (firstFrame == 1) && (lastFrame == -1)
    return % do nothing if not required
end

[~, ~, numFrames] = size(obj.M0_ff_video);

if firstFrame > 0 && firstFrame < numFrames || lastFrame > 1 && lastFrame <= numFrames

    tic

    fprintf("    - Video Cropping from %d to %d frames...\n", firstFrame, lastFrame);

    if lastFrame == -1
        lastFrame = numFrames;
    end

    if firstFrame == -1
        firstFrame = 1;
    end

    obj.M0_ff_video = obj.M0_ff_video(:, :, firstFrame:lastFrame);
    obj.M0_data_video = obj.M0_data_video(:, :, firstFrame:lastFrame);
    obj.M1_data_video = obj.M1_data_video(:, :, firstFrame:lastFrame);
    obj.M2_data_video = obj.M2_data_video(:, :, firstFrame:lastFrame);

    if ~isempty(obj.SH_data_hypervideo)
        obj.SH_data_hypervideo = obj.SH_data_hypervideo(:, :, :, firstFrame:lastFrame);
    end

    disp(['Data cube frame: ', num2str(firstFrame), '/', num2str(numFrames), ' to ', num2str(lastFrame), '/', num2str(numFrames)])

    fprintf("    - Video Cropping took: %ds\n", round(toc));

else
    disp('Wrong value for the first frame. Set as 1.')
    disp('Wrong value for the last frame. Set as the end.')
    disp(['Data cube frame: 1/', num2str(numFrames), ' to ', num2str(numFrames), '/', num2str(numFrames)])
end

end
