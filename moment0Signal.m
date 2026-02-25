function moment0Signal(M0)
ToolBox = getGlobalToolBox;
params = ToolBox.getParams;
jsonParams = params.json;


maskArtery = ToolBox.Cache.maskArtery;
maskVein = ToolBox.Cache.maskVein;
maskVessel = ToolBox.Cache.maskVessel;
xy_barycenter = ToolBox.Cache.xy_barycenter;

% Validating inputs
if ~any(maskArtery)
    error("Given Mask Artery is empty.")
end

[numX, numY, numFrames] = size(M0);
% Create section mask
r1 = jsonParams.SizeOfField.SmallRadiusRatio;
r2 = jsonParams.SizeOfField.BigRadiusRatio;
maskSection = diskMask(numX, numY, r1, r2);

% Create vessel masks
maskArterySection = imresize(maskArtery,[numX, numY]) & maskSection;


signal = squeeze(sum(M0.*maskArterySection,[1,2])/nnz(maskArterySection));

ToolBox.Output.add("ArteryM0Signal",signal,h5path = "Artery/Velocity/m0SectionSignalRaw/Signal");
ToolBox.Output.add("ArteryM0SignalnnzSection",nnz(maskArterySection),h5path = "Artery/Velocity/m0SectionSignalRaw/nnzSection");
end