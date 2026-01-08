function VelocityStatisticalOutputs(v_RMS_video, maskArtery, maskVein, maskArterySection, maskVeinSection)

% Initial Setup
ToolBox = getGlobalToolBox;
params = ToolBox.getParams;

mask = NaN(size(maskArterySection));
mask(maskArterySection) = 1;
maskArterySectionNaN = mask;

mask = NaN(size(maskVeinSection));
mask(maskVeinSection) = 1;
maskVeinSectionNaN = mask;

ArteryMeanSection = mean(v_RMS_video .* (maskArterySectionNaN), 'all', 'omitnan');
ArteryMaxSection = max(v_RMS_video .* (maskArterySectionNaN), [], 'all', 'omitnan');
ArteryMinSection = min(v_RMS_video .* (maskArterySectionNaN), [], 'all', 'omitnan');

ToolBox.Output.add('ArteryVelocityMeanSection', ArteryMeanSection, 'mm/s', h5path = '/Artery/Velocity/MeanSection');
ToolBox.Output.add('ArteryVelocityMaxSection', ArteryMaxSection, 'mm/s', h5path = '/Artery/Velocity/MaxSection');
ToolBox.Output.add('ArteryVelocityMinSection', ArteryMinSection, 'mm/s', h5path = '/Artery/Velocity/MinSection');

VeinMeanSection = mean(v_RMS_video .* (maskVeinSectionNaN), 'all', 'omitnan');
VeinMaxSection = max(v_RMS_video .* (maskVeinSectionNaN), [], 'all', 'omitnan');
VeinMinSection = min(v_RMS_video .* (maskVeinSectionNaN), [], 'all', 'omitnan');

ToolBox.Output.add('VeinVelocityMeanSection', VeinMeanSection, 'mm/s', h5path = '/Vein/Velocity/MeanSection');
ToolBox.Output.add('VeinVelocityMaxSection', VeinMaxSection, 'mm/s', h5path = '/Vein/Velocity/MaxSection');
ToolBox.Output.add('VeinVelocityMinSection', VeinMinSection, 'mm/s', h5path = '/Vein/Velocity/MinSection');

% Mode and median for artery section
ArteryModeSection = mode(v_RMS_video .* (maskArterySectionNaN), 'all'); % omits nan by default
ArteryMedianSection = median(v_RMS_video .* (maskArterySectionNaN), 'all', 'omitnan');

ToolBox.Output.add('ArteryVelocityModeSection', ArteryModeSection, 'mm/s', h5path = '/Artery/Velocity/ModeSection');
ToolBox.Output.add('ArteryVelocityMedianSection', ArteryMedianSection, 'mm/s');

% Mode and median for vein section
VeinModeSection = mode(v_RMS_video .* maskVeinSectionNaN, 'all'); % omits nan by default
VeinMedianSection = median(v_RMS_video .* maskVeinSectionNaN, 'all', 'omitnan');

ToolBox.Output.add('VeinVelocityModeSection', VeinModeSection, 'mm/s', h5path = '/Vein/Velocity/ModeSection');
ToolBox.Output.add('VeinVelocityMedianSection', VeinMedianSection, 'mm/s', h5path = '/Vein/Velocity/MedianSection');

% Remove 25-75 percentiles for artery section
arteryData = v_RMS_video .* (maskArterySectionNaN);
arteryData(isnan(arteryData)) = []; % remove all NaN

q25_artery = prctile(arteryData, 25);
q75_artery = prctile(arteryData, 75);
ToolBox.Output.add('ArteryVelocity25Percentile', q25_artery, 'mm/s', h5path = '/Artery/Velocity/Percentile25');
ToolBox.Output.add('ArteryVelocity75Percentile', q75_artery, 'mm/s', h5path = '/Artery/Velocity/Percentile75');

arteryTrimmed = arteryData(arteryData >= q25_artery & arteryData <= q75_artery);

ArteryMeanTrimmed = mean(arteryTrimmed); % already no nans
ArteryMedianTrimmed = median(arteryTrimmed);
ArteryModeTrimmed = mode(arteryTrimmed);

ToolBox.Output.add('ArteryVelocityMeanTrimmed', ArteryMeanTrimmed, 'mm/s', h5path = '/Artery/Velocity/MeanTrimmed');
ToolBox.Output.add('ArteryVelocityMedianTrimmed', ArteryMedianTrimmed, 'mm/s', h5path = '/Artery/Velocity/MedianTrimmed');
ToolBox.Output.add('ArteryVelocityModeTrimmed', ArteryModeTrimmed, 'mm/s', h5path = '/Artery/Velocity/ModeTrimmed');

% Remove 25-75 percentiles for vein section
veinData = v_RMS_video .* maskVeinSectionNaN;
veinData(isnan(veinData)) = []; % remove all NaN

q25_vein = prctile(veinData, 25);
q75_vein = prctile(veinData, 75);
ToolBox.Output.add('VeinVelocity25Percentile', q25_vein, 'mm/s', h5path = '/Vein/Velocity/Percentile25');
ToolBox.Output.add('VeinVelocity75Percentile', q75_vein, 'mm/s', h5path = '/Vein/Velocity/Percentile75');
veinTrimmed = veinData(veinData >= q25_vein & veinData <= q75_vein);

VeinMeanTrimmed = mean(veinTrimmed, 'omitnan');
VeinMedianTrimmed = median(veinTrimmed);
VeinModeTrimmed = mode(veinTrimmed);

ToolBox.Output.add('VeinVelocityMeanTrimmed', VeinMeanTrimmed, 'mm/s', h5path = '/Vein/Velocity/MeanTrimmed');
ToolBox.Output.add('VeinVelocityMedianTrimmed', VeinMedianTrimmed, 'mm/s', h5path = '/Vein/Velocity/MedianTrimmed');
ToolBox.Output.add('VeinVelocityModeTrimmed', VeinModeTrimmed, 'mm/s', h5path = '/Vein/Velocity/ModeTrimmed');

end
