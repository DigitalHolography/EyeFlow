classdef Cache < handle
% Cache class used to cache small variables through the execution
% Used for the handle
properties
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
    fs double % cached frame rate
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
    HeartBeatFFT double % cached heartbeat frequency
    HeartBeatFFTSTE double % cached heartbeat frequency standard error
    harmonics double % cached harmonics frequencies

end

methods

    function obj = Cache()
        obj.createColorMaps();
    end

    function createColorMaps(obj)
        obj.cmapArtery = cmapLAB(256, [0 0 0], 0, [1 0 0], 1/3, [1 1 0], 2/3, [1 1 1], 1);
        obj.cmapVein = cmapLAB(256, [0 0 0], 0, [0 0 1], 1/3, [0 1 1], 2/3, [1 1 1], 1);
        obj.cmapAV = cmapLAB(256, [0 0 0], 0, [1 0 1], 1/3, [1 1 1], 1);
    end

end

end
