classdef CacheClass < handle
% Cache class used to cache small variables through the execution
% Used for the handle
properties

    % Time Vector
    t double %cached time vector

    % Masks
    maskArtery logical % cached mask artery
    maskVein logical % cached mask vein
    maskNeighbors logical % cached mask neighbors
    maskBackground logical % cached mask background
    scoreMaskArtery double % cached score mask artery
    scoreMaskVein double % cached score mask vein

    % Diastolic and systolic indices
    sysIdxList double % cached systolic indices list
    sysIdx double % cached systolic index
    diasIdx double % cached diastolic index

    % Other cached variables
    papillaDiameter double % cached papilla diameter
    M0_RGB double % cached mean image

    % Color maps
    cmapArtery double % cached colormap artery
    cmapVein double % cached colormap vein
    cmapAV double % cached colormap AV

    % Barycenter
    xy_barycenter double % cached xy_barycenter
    xy_CRA double % cached xy_CRA
    xy_CRV double % cached xy_CRV
    xy_papilla double

    % Heartbeat
    HeartBeatFFT double % cached heartbeat frequency in Hz
    HeartBeatFFTSTE double % cached heartbeat frequency standard error in Hz
    harmonics double % cached harmonics frequencies

end

methods

    function obj = CacheClass()
        obj.createColorMaps();
    end

    function createColorMaps(obj)
        obj.cmapArtery = cmapLAB(256, [0 0 0], 0, [1 0 0], 1/3, [1 1 0], 2/3, [1 1 1], 1);
        obj.cmapVein = cmapLAB(256, [0 0 0], 0, [0 0 1], 1/3, [0 1 1], 2/3, [1 1 1], 1);
        obj.cmapAV = cmapLAB(256, [0 0 0], 0, [1/2 0 1/2], 1/3, [1 0 1], 2/3, [1 1 1], 1);
    end

    function createtimeVector(obj, ToolBox, numFrames)
        time_stamps = ToolBox.record_time_stamps_us;
        params = ToolBox.getParams;
        startFrame = params.json.Preprocess.Crop.StartFrame;
        endFrame = params.json.Preprocess.Crop.EndFrame;

        if ~isempty(time_stamps) ...
                && ~isempty(ToolBox.holo_frames) ...
                && params.json.use_time_stamps
            % if record_time_stamps_us
            % was specified in the HD folder from the .holo footer

            % binsize = 64;
            % (ToolBox.holo_frames.last - ToolBox.holo_frames.first + 1) / ToolBox.stride;
            if time_stamps.last <= time_stamps.first
                fprintf('Invalid time stamps detected, using default fs and stride');
                time_array = getTimeLinear(ToolBox, numFrames);
            else

                try
                    % Create time vector from time stamps
                    time_array = getTimeTimestamp(time_stamps, numFrames);
                catch
                    % If error occurs, use default fs and stride
                    fprintf('Could not create time vector from time stamps, using default fs and stride');
                    time_array = getTimeLinear(ToolBox, numFrames);
                end

            end

        else
            % Create time vector from fs and stride
            time_array = getTimeLinear(ToolBox, numFrames);
        end

        if startFrame < 1
            startFrame = 1;
        end

        if endFrame > numFrames || endFrame == -1
            endFrame = numFrames;
        end

        obj.t = time_array(startFrame:endFrame);

    end

end

end

function time_array = getTimeTimestamp(time_stamps, numFrames)
t1 = 0;
t2 = (time_stamps.last - time_stamps.first) / 1e6; % in s
time_array = linspace(t1, t2, numFrames);
end

function time_array = getTimeLinear(ToolBox, numFrames)
t1 = 0;
t2 = ToolBox.stride / (ToolBox.fs * 1000) * (numFrames - 1); % in s
time_array = linspace(t1, t2, numFrames);
end
