classdef OutputClass < handle
% Class to hold the Output of the retinal flow analysis pipeline

properties
    % % General
    % NumFrames
    % FrameRate
    % InterFramePeriod

    % % Arterial waveform analysis
    % SystoleIndices
    % HeartBeat
    % ArterySystoleMaxIndices
    % ArteryDiastoleMinIndices
    % ArteryTimeToMaxIncrease
    % ArteryTimeToPeakSystole
    % ArteryTimeToMinDiastole
    % ArteryTimeToPeakSystoleFromDiastole
    % ArteryTimeToDescent
    % ArteryTimePeakToDescent

    % % Venous waveform analysis
    % VeinTimeToPeakFromMin

    % % Velocity (mm/s)
    % ArteryVelocityMean
    % ArteryVelocityMax
    % ArteryVelocityMin
    % VeinVelocityMean
    % VeinVelocityMax
    % VeinVelocityMin

    % % Velocity Statistical Outputs (mm/s)
    % ArteryVelocityMeanSection
    % ArteryVelocityMaxSection
    % ArteryVelocityMinSection
    % ArteryVelocityModeSection
    % ArteryVelocityMedianSection
    % ArteryVelocityMeanTrimmed
    % ArteryVelocityMedianTrimmed
    % ArteryVelocityModeTrimmed

    % VeinVelocityMeanSection
    % VeinVelocityMaxSection
    % VeinVelocityMinSection
    % VeinVelocityModeSection
    % VeinVelocityMedianSection
    % VeinVelocityMeanTrimmed
    % VeinVelocityMedianTrimmed
    % VeinVelocityModeTrimmed

    % ArteryVelocity25Percentile
    % ArteryVelocity75Percentile
    % VeinVelocity25Percentile
    % VeinVelocity75Percentile

    % % Flow Rate (µL/s or mm³/s)
    % ArteryFlowRateMean
    % ArteryFlowRateMax
    % ArteryFlowRateMin
    % VeinFlowRateMean
    % VeinFlowRateMax
    % VeinFlowRateMin

    % % Resistivity and Pulsatility Indices (dimensionless)
    % ArteryResistivityIndexVelocity
    % ArteryPulsatilityIndexVelocity
    % ArteryMaxMinRatioVelocity
    % VeinResistivityIndexVelocity
    % VeinPulsatilityIndexVelocity
    % VeinMaxMinRatioVelocity

    % ArteryResistivityIndexFlowRate
    % ArteryPulsatilityIndexFlowRate
    % ArteryMaxMinRatioFlowRate
    % VeinResistivityIndexFlowRate
    % VeinPulsatilityIndexFlowRate
    % VeinMaxMinRatioFlowRate

    % % Stroke Volume & Cycle fractions
    % ArteryCycleVolume
    % ArterySystolicFraction
    % ArteryDiastolicFraction
    % VeinCycleVolume
    % VeinSystolicFraction
    % VeinDiastolicFraction

    % % Temporal durations
    % SystoleDuration
    % DiastoleDuration
    % SystolicUpstroke
    % SystolicDownstroke
    % DiastolicRunoff

    % % Vessel Diameters
    % ArteryDiameterMean
    % ArteryDiameterMedian
    % ArteryDiameterSpread
    % ArteryValidSectionCount
    % ArteryTotalSectionCount
    % VeinDiameterMean
    % VeinDiameterMedian
    % VeinDiameterSpread
    % VeinValidSectionCount
    % VeinTotalSectionCount

    % % Additional hemodynamic parameters
    % ArteryPulseWaveVelocity
    % VeinPulseWaveVelocity
    % DicroticNotchVisibility
    % ViscosityDuringSystole
    % ViscosityDuringDiastole
    % ViscosityMean
    % PapillaRatio
    % WindkesselDecayRC
    % WindkesselPureDelay
    % WindkesselCeff
    % WindkesselReff
    % ArteryVeinPhaseDelay

    % % Mask
    % ArteryArea % in mm²
    % VeinArea
    % RemainingArea
    % ArteryNbPxl % in pixels on all field
    % VeinNbPxl
    % RemainingNbPxl
    % ArterySelNbPxl % in pixels on targeted region (smallRadius bigRadius)
    % VeinSelNbPxl
    % RemainingSelNbPxl

    % % Mic
    % PredictedEyeSide

    % % Quality Control Scores
    % QualityControlScoreMaskArtery
    % QualityControlScoreMaskVein

    % % Time info
    % UnixTimestampFirst
    % UnixTimestampLast

    data
end

methods

    function obj = OutputClass()
    end

    function add(obj, name, value, unit, ste, vars)

        arguments
            obj
            name string
            value
            unit = "" % old behavior
            ste = NaN % old behavior
            vars.h5path string = ""
            vars.unit string = ""
            vars.standard_error = NaN
        end

        if nargin > 3
            vars.unit = unit;
        end

        if nargin > 4
            vars.standard_error = ste;
        end

        obj.data.(name).value = value;
        obj.data.(name).standard_error = vars.standard_error;
        obj.data.(name).unit = vars.unit;

        if vars.h5path == ""
            obj.data.(name).h5path = sprintf("/%s", name);
        else
            obj.data.(name).h5path = (vars.h5path);
        end

    end

    function writeJson(obj, path)
        props = fieldnames(obj.data);

        d = containers.Map('KeyType', 'char', 'ValueType', 'any');

        for i = 1:length(props)
            v = obj.data.(props{i}).value;

            if isscalar(v)

                if obj.data.(props{i}).unit == ""
                    key = props{i};
                else
                    key = props{i} + "_" + obj.data.(props{i}).unit;
                end

                d(char(key)) = v;
            end

        end

        if ~isempty(d)
            jsonText = jsonencode(d, "PrettyPrint", true);

            fid = fopen(path, 'w');

            if fid == -1
                error('Cannot open file for writing: %s', path);
            end

            fwrite(fid, jsonText, 'char');
            fclose(fid);
        end

    end

    function writeHdf5(obj, path)

        [folder_dir, folder_name, ~] = fileparts(path);
        file_path = fullfile(folder_dir, strcat(folder_name, ".h5"));

        if isfile(file_path)
            delete(file_path)
        end

        props = fieldnames(obj.data);

        for i = 1:length(props)

            h5path = (obj.data.(props{i}).h5path);
            temp = char(h5path);

            if temp(1) ~= '/'
                h5path = strcat("/", h5path);
            end

            if ~isnan(obj.data.(props{i}).standard_error)
                writeNumericToHDF5(file_path, strcat(h5path, "/value"), obj.data.(props{i}).value);
                h5writeatt(file_path, strcat(h5path, "/value"), "unit", obj.data.(props{i}).unit);
                writeNumericToHDF5(file_path, strcat(h5path, "/ste"), obj.data.(props{i}).standard_error);
                h5writeatt(file_path, strcat(h5path, "/ste"), "unit", obj.data.(props{i}).unit);
            else
                writeNumericToHDF5(file_path, h5path, obj.data.(props{i}).value);
                h5writeatt(file_path, h5path, "unit", obj.data.(props{i}).unit);
            end

            h5writeatt(file_path, h5path, "nameID", props{i});

        end

    end

end

end

% --- Helper: write numeric dataset ---

function writeNumericToHDF5(path, datasetPath, value)

if ~isempty(value)

    if isnumeric(value) & ~isreal(value)
        warning("Complex values should be handled before call : %s", datasetPath);
        return;
    end

    if islogical(value) % Hdf5 matlab does not support
        value = uint8(value);
    end

    stored_value = whos('value');
    bytesize = stored_value.bytes;
    original_class = stored_value.class;

    % Threshold for precision reduction
    threshold = 1e6; % 1 MBytes

    % Reduce precision if numeric type allows it
    isReduced = false;

    if bytesize > threshold
        isReduced = true;
        mini = min(value(:));
        MM = max(value(:));
        value = uint8(rescale(value) * 255);
    end

    h5create(path, datasetPath, size(value), 'Datatype', class(value), 'Chunksize', size(value), 'Deflate', 6);
    h5write(path, datasetPath, value);

    if isReduced
        h5writeatt(path, datasetPath, "minimum", mini);
        h5writeatt(path, datasetPath, "maximum", MM);
        h5writeatt(path, datasetPath, "original_class", original_class);
    end

end

end
