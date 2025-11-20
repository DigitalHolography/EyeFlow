classdef AINetworksClass < handle
% Handle class to store AI models net one and only once
properties
    VesselSegmentationParam
    VesselSegmentationNet
    AVSegmentationParam
    AVSegmentationNet
    OpticDiskDetectorParam
    OpticDiskDetectorNet
end

methods

    function AINetworksObj = AINetworksClass()
        % Constructor initializes empty properties
        AINetworksObj.VesselSegmentationParam = '';
        AINetworksObj.VesselSegmentationNet = [];
        AINetworksObj.AVSegmentationParam = '';
        AINetworksObj.AVSegmentationNet = [];
        AINetworksObj.OpticDiskDetectorParam = false;
        AINetworksObj.OpticDiskDetectorNet = [];
    end

    function updateAINetworks(obj, params)
        % Load or update AI networks based on parameters
        % Check for change to avoid redundant loading

        if ~isfolder('Models')
            mkdir('Models');
        end

        maskParams = params.json.Mask;

        if ~strcmp(maskParams.VesselSegmentationMethod, 'AI') && ...
                ~maskParams.AVDiasysSegmentationNet && ...
                ~maskParams.AVCorrelationSegmentationNet && ...
                ~maskParams.OpticDiskDetectorNet

            fprintf("    - No AI Networks selected in parameters. Skipping loading.\n");
            return; % No AI networks needed
        end

        % Load Optic Disk Detector Network if parameter changed
        OpticDiskDetectorParamChanged = ~isequal(obj.OpticDiskDetectorParam, maskParams.OpticDiskDetectorNet);

        if OpticDiskDetectorParamChanged && maskParams.OpticDiskDetectorNet
            tic
            fprintf("    - Loading Optic Disk Detector Network...\n");
            [obj.OpticDiskDetectorNet] = loadOpticDiskNetwork();
            obj.OpticDiskDetectorParam = maskParams.OpticDiskDetectorNet;
            fprintf("    - Loading Optic Disk Detector Network took: %ds\n", round(toc));
        end

        % Load Vessel Segmentation Network if parameter changed
        VesselSegmentationParamChanged = ~isequal(obj.VesselSegmentationParam, maskParams.VesselSegmentationMethod);

        if VesselSegmentationParamChanged && strcmp(maskParams.VesselSegmentationMethod, 'AI')
            tic
            fprintf("    - Loading Vessel Segmentation Networks...\n");
            [obj.VesselSegmentationNet] = loadVesselSegmentationNetworks();
            obj.VesselSegmentationParam = maskParams.VesselSegmentationMethod;
            fprintf("    - Loading Vessel Segmentation Networks took: %ds\n", round(toc));
        end

        % Determine AV Segmentation Parameter
        if maskParams.AVDiasysSegmentationNet && maskParams.AVCorrelationSegmentationNet
            AVSegParam = 'AVBothSegmentationNet';
        elseif maskParams.AVDiasysSegmentationNet
            AVSegParam = 'AVDiasysSegmentationNet';
        elseif maskParams.AVCorrelationSegmentationNet
            AVSegParam = 'AVCorrelationSegmentationNet';
        else
            AVSegParam = '';
        end

        % Load AV Segmentation Network if parameter changed
        AVSegmentationParamChanged = ~isequal(obj.AVSegmentationParam, AVSegParam);

        if AVSegmentationParamChanged
            tic
            fprintf("    - Loading AV Segmentation Networks...\n");
            [obj.AVSegmentationNet] = loadAVSegmentationNetworks(params);
            obj.AVSegmentationParam = AVSegParam;
            fprintf("    - Loading AV Segmentation Networks took: %ds\n", round(toc));
        end

        if ~AVSegmentationParamChanged && ~VesselSegmentationParamChanged
            fprintf("    - AI Networks are up to date. No loading needed.\n");
        end

    end

end

end
