function [VesselnessNet, VesselNet] = loadAINetworks(params)
% loadAINetworks Loads AI networks for vesselness and artery/vein segmentation
%   Inputs:
%       params - Parameters structure containing segmentation settings
%   Outputs:
%       VesselnessNet - Loaded vesselness segmentation network
%       VesselNet - Loaded artery/vein segmentation network

if ~isfolder('Models')
    mkdir('Models');
end

VesselnessNet = [];
VesselNet = [];

if strcmp(params.json.Mask.VesselSegmentationMethod, 'AI')

    model_name = "iternet5_vesselness";
    model_path = 'Models\iternet5_vesselness.onnx';

    if ~isfile(model_path)
        % Download the model from Hugging Face
        url = sprintf('https://huggingface.co/DigitalHolography/%s/resolve/main/%s', model_name, model_name);
        fprintf(url)
        websave(model_path, url);
    end

    % Import the ONNX network
    try
        % Try the newer function first
        VesselnessNet = importNetworkFromONNX(model_path);
    catch
        % Fall back to the older function
        warning('off')
        VesselnessNet = importONNXNetwork(model_path);
        warning('on')
    end

    fprintf("Loaded iternet5_vesselness model for vessel segmentation\n");

end

if params.json.Mask.AVDiasysSegmentationNet

    model_name = "iternet5_av_diasys";
    model_path = 'Models\iternet5_av_diasys.onnx';

    if ~isfile(model_path)
        % Download the model from Hugging Face
        url = sprintf('https://huggingface.co/DigitalHolography/%s/resolve/main/%s', model_name);
        websave(model_path, url);
    end

    try
        % Try the newer function first
        VesselNet = importNetworkFromONNX('Models\iternet5_av_diasys.onnx');
    catch
        % Fall back to the older function
        warning('off')
        VesselNet = importONNXNetwork('Models\iternet5_av_diasys.onnx');
        warning('on')
    end

    fprintf("Loaded iternet5_av_diasys model for artery/vein segmentation\n");

elseif params.json.Mask.AVCorrelationSegmentationNet

    try
        % Try the newer function first
        VesselNet = importNetworkFromONNX('Models\iternet_5_av_corr.onnx');
    catch
        % Fall back to the older function
        warning('off')
        VesselNet = importONNXNetwork('Models\iternet_5_av_corr.onnx');
        warning('on')
    end

    fprintf("Loaded iternet_5_av_corr model for artery/vein segmentation\n");

elseif params.json.Mask.AVCorrelationSegmentationNet && params.json.Mask.AVDiasysSegmentationNet

    try
        % Try the newer function first
        VesselNet = importNetworkFromONNX('Models\iternet5_av_corr_diasys.onnx');
    catch
        % Fall back to the older function
        warning('off')
        VesselNet = importONNXNetwork('Models\iternet5_av_corr_diasys.onnx');
        warning('on')
    end

    fprintf("Loaded iternet5_av_corr_diasys model for artery/vein segmentation\n");

end

end
