classdef OutputClass < handle
% Class to hold the Output of the retinal flow analysis pipeline

properties
    % General
    NumFrames
    FrameRate
    InterFramePeriod

    % Arterial waveform analysis
    SystoleIndices
    HeartBeat
    ArterySystoleMaxIndices
    ArteryDiastoleMinIndices
    ArteryTimeToMaxIncrease
    ArteryTimeToPeakSystole
    ArteryTimeToMinDiastole
    ArteryTimeToPeakSystoleFromDiastole
    ArteryTimeToDescent
    ArteryTimePeakToDescent

    % Venous waveform analysis
    VeinTimeToPeakFromMin

    % Velocity (mm/s)
    ArteryVelocityMean
    ArteryVelocityMax
    ArteryVelocityMin
    VeinVelocityMean
    VeinVelocityMax
    VeinVelocityMin

    % Flow Rate (µL/s or mm³/s)
    ArteryFlowRateMean
    ArteryFlowRateMax
    ArteryFlowRateMin
    VeinFlowRateMean
    VeinFlowRateMax
    VeinFlowRateMin

    % Resistivity and Pulsatility Indices (dimensionless)
    ArteryResistivityIndexVelocity
    ArteryPulsatilityIndexVelocity
    ArteryMaxMinRatioVelocity
    VeinResistivityIndexVelocity
    VeinPulsatilityIndexVelocity
    VeinMaxMinRatioVelocity

    ArteryResistivityIndexFlowRate
    ArteryPulsatilityIndexFlowRate
    ArteryMaxMinRatioFlowRate
    VeinResistivityIndexFlowRate
    VeinPulsatilityIndexFlowRate
    VeinMaxMinRatioFlowRate

    % Stroke Volume & Cycle fractions
    ArteryCycleVolume
    ArterySystolicFraction
    ArteryDiastolicFraction
    VeinCycleVolume
    VeinSystolicFraction
    VeinDiastolicFraction

    % Temporal durations
    SystoleDuration
    DiastoleDuration
    SystolicUpstroke
    SystolicDownstroke
    DiastolicRunoff

    % Vessel Diameters
    ArteryDiameterMean
    ArteryDiameterMedian
    ArteryDiameterSpread
    ArteryValidSectionCount
    ArteryTotalSectionCount
    VeinDiameterMean
    VeinDiameterMedian
    VeinDiameterSpread
    VeinValidSectionCount
    VeinTotalSectionCount

    % Additional hemodynamic parameters
    ArteryPulseWaveVelocity
    VeinPulseWaveVelocity
    DicroticNotchVisibility
    ViscosityDuringSystole
    ViscosityDuringDiastole
    ViscosityMean
    PapillaRatio
    WindkesselDecayRC
    WindkesselPureDelay
    WindkesselCeff
    WindkesselReff
    ArteryVeinPhaseDelay

    % Time info
    UnixTimestampFirst
    UnixTimestampLast

    % Signals & Extra & DimOut
    Signals SignalsClass
    Extra ExtraClass
    DimOut DimOutClass
end

methods

    function obj = OutputClass()
        % Constructor for the class, fills the properties with default values
        props = setdiff(properties(obj), "Signals");
        props = setdiff(props, "Extra");
        props = setdiff(props, "DimOut");

        for i = 1:length(props)
            obj.(props{i}).value = NaN;
            obj.(props{i}).standard_error = NaN;
            obj.(props{i}).unit = "";
        end

        obj.Signals = SignalsClass();
        obj.Extra = ExtraClass();
        obj.DimOut = DimOutClass();
    end

    function add(obj, name, value, unit, standard_error)

        if nargin < 5
            standard_error = NaN;
        end

        if isprop(obj, name)
            obj.(name).value = value;
            obj.(name).standard_error = standard_error;
            obj.(name).unit = unit;
        else
            error('Property %s does not exist in OutputClass', name);
        end

    end

    function writeJson(obj, path)
        props = setdiff(properties(obj), "Signals");
        props = setdiff(props, "Extra");
        props = setdiff(props, "DimOut");
        
        data = struct();

        for i = 1:length(props)
            data.(props{i}) = obj.(props{i}).value;
            data.(strcat(props{i}, "_ste")) = obj.(props{i}).standard_error;
            data.(strcat(props{i}, "_unit")) = obj.(props{i}).unit;
        end

        jsonText = jsonencode(data, "PrettyPrint", true);
        fid = fopen(path, 'w');

        if fid == -1
            error('Cannot open file for writing: %s', path);
        end

        fwrite(fid, jsonText, 'char');
        fclose(fid);
    end

    function writeHdf5(obj, path)
        props = setdiff(properties(obj), "Signals");
        props = setdiff(props, "Extra");
        props = setdiff(props, "DimOut");
        
        [folder_dir, folder_name, ~] = fileparts(path);
        file_path = fullfile(folder_dir, strcat(folder_name, ".h5"));

        if isfile(file_path)
            delete(file_path)
        end

        for i = 1:length(props)
            if ~isempty(obj.(props{i}).value)
                writeNumericToHDF5(file_path, strcat("/", "Scalars", "/", props{i}, "/", props{i}), obj.(props{i}).value);
                h5writeatt(file_path, strcat("/", "Scalars", "/", props{i}, "/", props{i}), "unit", obj.(props{i}).unit);
            end
            
            if ~isempty(obj.(props{i}).standard_error)
                writeNumericToHDF5(file_path, strcat("/", "Scalars", "/", props{i}, "/", props{i}, "_ste"), obj.(props{i}).standard_error);
                h5writeatt(file_path, strcat("/", "Scalars", "/", props{i}, "/", props{i}, "_ste"), "unit", obj.(props{i}).unit);
            end
            
            
        end

        % Save Signals and Extra in a separate directory of the same h5 file
        obj.Signals.writeHdf5(file_path);
        obj.Extra.writeHdf5(file_path);
        obj.DimOut.writeHdf5(file_path);

    end

end

end
