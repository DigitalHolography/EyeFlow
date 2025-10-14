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
        obj.cmapAV = cmapLAB(256, [0 0 0], 0, [1 0 1], 1/3, [1 1 1], 1);
    end

    function createtimeVector(obj, ToolBox, numFrames)
        time_stamps = ToolBox.record_time_stamps_us;
        t1 = 0;
        t2 = (time_stamps.last - time_stamps.first) / 1e6;
        t_stamp = linspace(t1, t2, numFrames - 1);
        t_stamp(end + 1) = t2 + t_stamp(2);
        obj.t = t_stamp;
    end

end

end
