classdef AINetworksClass < handle
% Handle class to store AI models net one and only once
properties
    VesselSegmentationParam
    VesselSegmentationNet
    AVSegmentationParam
    AVSegmentationNet
    OpticDiskDetectorParam
    OpticDiskDetectorNet
    EyeSideClassifierParam
    EyeSideClassifierNet
    EyeDiaphragmSegmentationParam
    EyeDiaphragmSegmentationNet
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

        pythonFolder = 'C:\Users\Michael\Documents\MATLAB\EyeFlow\Preprocessing';
        if count(py.sys.path, pythonFolder) == 0
            insert(py.sys.path, int32(0), pythonFolder);
        end

        py.importlib.invalidate_caches(); % refresh import cache

        maskParams = params.json.Mask;

        if ~strcmp(maskParams.VesselSegmentationMethod, 'AI') && ...
                ~maskParams.AVDiasysSegmentationNet && ...
                ~maskParams.AVCorrelationSegmentationNet && ...
                ~maskParams.OpticDiskDetectorNet && ...
                ~maskParams.PapillaDiskDetectorNet && ...
                ~maskParams.EyeSideClassifierNet

            fprintf("    - No AI Networks selected in parameters. Skipping loading.\n");
            return; % No AI networks needed
        end

        % Load Optic Disk Detector Network if parameter changed
        OpticDiskDetectorParamChanged = ~isequal(obj.OpticDiskDetectorParam, maskParams.OpticDiskDetectorNet);

        if OpticDiskDetectorParamChanged && maskParams.OpticDiskDetectorNet
            pyenv; % Ensure Python environment is initialized

            tic
            fprintf("    - Loading Optic Disk Detector Network...\n");
            % [obj.OpticDiskDetectorNet] = loadOpticDiskNetwork();
            [obj.OpticDiskDetectorNet] = py.import_yolo_model.import_yolo_model("runs/train/optic_disk_segmentation2/weights/best.pt"); % Load the trained model
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

        % Load Eye Side Classifier Network if parameter changed
        EyeSideClassifierParamChanged = ~isequal(obj.EyeSideClassifierParam, maskParams.EyeSideClassifierNet);

        if EyeSideClassifierParamChanged && maskParams.EyeSideClassifierNet
            tic
            fprintf("    - Loading Eye Side Classifier Network...\n");
            [obj.EyeSideClassifierNet] = loadEyeSideClassifierNetwork();
            obj.EyeSideClassifierParam = maskParams.EyeSideClassifierNet;
            fprintf("    - Loading Eye Side Classifier Network took: %ds\n", round(toc));
        end

        % Load Eye Diaphragm Segmentation Network if parameter changed
        EyeDiaphragmSegmentationParamChanged = ~isequal(obj.EyeDiaphragmSegmentationParam, maskParams.EyeDiaphragmSegmentationNet);

        if EyeDiaphragmSegmentationParamChanged && maskParams.EyeDiaphragmSegmentationNet
            tic
            fprintf("    - Loading Eye Diaphragm Segmentation Network...\n");
            [obj.EyeDiaphragmSegmentationNet] = loadEyeDiaphragmSegmentation();
            obj.EyeDiaphragmSegmentationParam = maskParams.EyeDiaphragmSegmentationNet;
            fprintf("    - Loading Eye Diaphragm Segmentation Network took: %ds\n", round(toc));
        end

        if ~AVSegmentationParamChanged ...
                && ~VesselSegmentationParamChanged ...
                && ~OpticDiskDetectorParamChanged ...
                && ~EyeSideClassifierParamChanged ...
                && ~PapillaDiskDetectorParamChanged
            fprintf("    - AI Networks are up to date. No loading needed.\n");
        end

    end

end

end
