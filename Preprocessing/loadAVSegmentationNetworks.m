function [AVSegmentationNet] = loadAVSegmentationNetworks(params)
% loadAVSegmentationNetworks Loads AI networks for vesselness and artery/vein segmentation
%   Inputs:
%       params - Parameters structure containing segmentation settings
%   Outputs:
%       VesselSegmentationNet - Loaded vesselness segmentation network

if params.json.Mask.AVDiasysSegmentationNet

    model_name = "iternet5_av_diasys";
    model_path = 'Models\iternet5_av_diasys.onnx';

    if ~isfile(model_path)
        % Download the model from Hugging Face
        url = sprintf('https://huggingface.co/DigitalHolography/%s/resolve/main/%s', model_name, model_name);
        websave(model_path, url);
    end

    try
        % Try the newer function first
        AVSegmentationNet = importNetworkFromONNX(model_path);
    catch
        % Fall back to the older function
        warning('off')
        AVSegmentationNet = importONNXNetwork(model_path);
        warning('on')
    end

    fprintf("Loaded iternet5_av_diasys model for artery/vein segmentation\n");

elseif params.json.Mask.AVCorrelationSegmentationNet

    model_name = "iternet5_av_corr";
    model_path = 'Models\iternet5_av_corr.onnx';

    if ~isfile(model_path)
        % Download the model from Hugging Face
        url = sprintf('https://huggingface.co/DigitalHolography/%s/resolve/main/%s', model_name, model_name);
        websave(model_path, url);
    end

    try
        % Try the newer function first
        AVSegmentationNet = importNetworkFromONNX(model_path);
    catch
        % Fall back to the older function
        warning('off')
        AVSegmentationNet = importONNXNetwork(model_path);
        warning('on')
    end

    fprintf("Loaded iternet_5_av_corr model for artery/vein segmentation\n");

elseif params.json.Mask.AVCorrelationSegmentationNet && params.json.Mask.AVDiasysSegmentationNet

    model_name = "iternet5_av_corr_diasys";
    model_path = 'Models\iternet5_av_corr_diasys.onnx';

    if ~isfile(model_path)
        % Download the model from Hugging Face
        url = sprintf('https://huggingface.co/DigitalHolography/%s/resolve/main/%s', model_name, model_name);
        websave(model_path, url);
    end

    try
        % Try the newer function first
        AVSegmentationNet = importNetworkFromONNX(model_path);
    catch
        % Fall back to the older function
        warning('off')
        AVSegmentationNet = importONNXNetwork(model_path);
        warning('on')
    end

    fprintf("Loaded iternet5_av_corr_diasys model for artery/vein segmentation\n");

end

end
