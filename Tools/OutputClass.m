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
    ArteryVeinPhaseDelay

    % Time info
    UnixTimestampFirst
    UnixTimestampLast

    % Signals
    Signals SignalsClass
end

methods

    function obj = OutputClass()
        % Constructor for the class, fills the properties with default values
        props = setdiff(properties(obj), "Signals");

        for i = 1:length(props)
            obj.(props{i}).value = NaN;
            obj.(props{i}).standard_error = NaN;
            obj.(props{i}).unit = "";
        end

        obj.Signals = SignalsClass();
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
        [dir, name, ~] = fileparts(path);
        path = fullfile(dir, strcat(name, ".h5"));

        if isfile(path)
            delete(path)
        end

        for i = 1:length(props)
            h5create(path, strcat("/", props{i}), size(obj.(props{i}).value));
            h5write(path, strcat("/", props{i}), obj.(props{i}).value);
            h5create(path, strcat("/", props{i}, "_ste"), size(obj.(props{i}).standard_error));
            h5write(path, strcat("/", props{i}, "_ste"), obj.(props{i}).standard_error);
            h5create(path, strcat("/", props{i}, "_unit"), [1 1], Datatype = "string");
            h5write(path, strcat("/", props{i}, "_unit"), string(obj.(props{i}).unit));
        end

        % Save Signals in a separate directory of the same h5 file
        obj.Signals.writeHdf5(path);
    end

end

end
