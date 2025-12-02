classdef PreprocessorClass < handle
% Handles all preprocessing operations

properties
    is_preprocessed = false
    M0_ff double
    M0 double
    M1 double
    M2 double
    SH double
    f_RMS double
    f_AVG double
    directory char
    filenames char
    param_name char
    displacementField

    firstFrameIdx double
    lastFrameIdx double
end

methods

    function obj = PreprocessorClass(directory, filenames, param_name)
        obj.directory = directory;
        obj.filenames = filenames;
        obj.param_name = param_name;
    end

    function preprocess(obj, executionObj)

        if any(isnan(executionObj.M0), 'all')
            error('NaN values found in M0 data. Please check the input file.');
        end

        params = Parameters_json(obj.directory, obj.param_name);

        obj.M0 = executionObj.M0;
        obj.M1 = executionObj.M1;
        obj.M2 = executionObj.M2;
        obj.SH = executionObj.SH;

        % Execute preprocessing steps
        obj.register(params);
        obj.crop(params);
        obj.normalizeLocally(params);
        obj.resize(params);
        obj.nonRigidRegister(params);
        obj.interpolate(params);
        obj.removeOutliers(params);

        obj.is_preprocessed = true;
    end

end

methods (Access = private)

    function register(obj, params)
        firstFrame = params.json.Preprocess.Register.StartFrame;
        lastFrame = params.json.Preprocess.Register.EndFrame;
        enableRegistration = params.json.Preprocess.Register.Enable;

        if (firstFrame == 1) && (lastFrame == -1) || ~enableRegistration
            return % do nothing if not required
        end

        tic
        fprintf("    - Video Registering...\n");
        % Rigid registration implementation
        VideoRegistering(obj, firstFrame, lastFrame);
        fprintf("    - Video Registration took: %ds\n", round(toc));
    end

    function crop(obj, params)

        firstFrame = params.json.Preprocess.Crop.StartFrame;
        lastFrame = params.json.Preprocess.Crop.EndFrame;

        if (firstFrame == 1) && (lastFrame == -1)
            return % do nothing if not required
        end

        tic
        fprintf("    - Video Cropping...\n");
        % Cropping implementation
        [firstFrameOut, lastFrameOut] = VideoCropping(obj, firstFrame, lastFrame);

        if firstFrameOut ~= firstFrame
            obj.firstFrameIdx = firstFrameOut;
        end

        if lastFrameOut ~= lastFrame
            obj.lastFrameIdx = lastFrameOut;
        end

        fprintf("    - Video Cropping took: %ds\n", round(toc));
    end

    function normalizeLocally(obj, params)
        tic
        fprintf("    - Local Normalization...\n");
        % Local normalization implementation
        VideoNormalizingLocally(obj, params);
        fprintf("    - Local Normalization took: %ds\n", round(toc));
    end

    function resize(obj, params)
        tic
        fprintf("    - Video Resizing...\n");
        % Resizing implementation
        VideoResizing(obj, params);
        fprintf("    - Video Resizing took: %ds\n", round(toc));
    end

    function nonRigidRegister(obj, params)

        if ~params.json.Preprocess.NonRigidRegisteringFlag
            return
        end

        apply = params.json.Preprocess.NonRigidRegisteringApply;

        tic
        fprintf("    - Video Non-Rigid Registration...\n");
        % Non-rigid registration implementation
        VideoNonRigidRegistering(obj, apply);
        fprintf("    - Video Non-Rigid Registration took: %ds\n", round(toc));
    end

    function interpolate(obj, params)
        kInterp = params.json.Preprocess.InterpolationFactor;

        if kInterp == 0
            return
        end

        tic
        fprintf("    - Video Interpolating...\n");
        % Interpolation implementation
        VideoInterpolating(obj, params);
        fprintf("    - Video Interpolation took: %ds\n", round(toc));
    end

    function removeOutliers(obj, params)

        if ~params.json.Preprocess.RemoveOutliers.Flag
            return
        end

        tic
        fprintf("    - Video Outlier Removal...\n");
        % Outlier removal implementation
        VideoRemoveOutliers(obj, params);
        fprintf("    - Video Outlier Removal took: %ds\n", round(toc));
    end

end

end
