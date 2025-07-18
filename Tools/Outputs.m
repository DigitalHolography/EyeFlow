classdef Outputs < handle
% Class to hold the main outputs of the retinal flow analysis pipeline

properties

    NumFrames
    FrameRate
    InterFramePeriod

    % Arterial Wave form analysis
    SystoleIndices
    HeartBeat
    MaximumSystoleIndices
    MinimumDiastoleIndices
    TimeToMaxIncreaseSystolic
    TimeToPeakSystole
    TimeToMinimumDiastole
    TimeToPeakSystoleFromMinimumDiastole
    TimeToDescent
    TimePeakToDescent

    % Vein Wave form analysis
    TimetoPeakFromMinVein
    TimetoAscentFromMinVein
    TimetoDescentToMinVein

    % Phase Delay
    PhaseDelay

    % Velocity
    ArterialMeanVelocity
    ArterialMaximumVelocity
    ArterialMinimumVelocity
    VenousMeanVelocity
    VenousMaximumVelocity
    VenousMinimumVelocity

    % Volume Rate
    ArterialMeanVolumeRate
    ArterialMaximumVolumeRate
    ArterialMinimumVolumeRate
    VenousMeanVolumeRate
    VenousMaximumVolumeRate
    VenousMinimumVolumeRate

    % Resistivity and Pulsatility Indices
    ArterialResistivityIndexVelocity
    ArterialPulsatilityIndexVelocity
    ArterialMaxMinRatioVelocity
    VenousResistivityIndexVelocity
    VenousPulsatilityIndexVelocity
    VenousMaxMinRatioVelocity
    ArterialResistivityIndexVolumeRate
    ArterialPulsatilityIndexVolumeRate
    ArterialMaxMinRatioVolumeRate
    VenousResistivityIndexVolumeRate
    VenousPulsatilityIndexVolumeRate
    VenousMaxMinRatioVolumeRate

    % Stroke Volume
    ArterialCycleVolume
    ArterialSystolicFraction
    ArterialDiastolicFraction
    VenousCycleVolume
    VenousSystolicFraction
    VenousDiastolicFraction

    SystoleDuration
    DiastoleDuration
    SystolicUpstroke
    SystolicDownstroke
    DiastolicRunoff

    % Vessel Diameters
    ArterialDiameterAverage
    ArterialDiameterMedian
    ArterialDiameterSpread
    ArterialValidSections
    TotalArterialSections
    VenousDiameterAverage
    VenousDiameterMedian
    VenousDiameterSpread
    VenousValidSections
    TotalVenousSections

    % Extra
    PulseWaveVelocity
    % DicroticNotchVisibility
    DynamicViscosityDuringSystole
    DynamicViscosityDuringDiastole
    DynamicViscosityAverage
    PapillaRatio

end

methods

    function obj = Outputs()
    end

    function initOutputs(obj)
        % Constructor for the class, fills the properties with default values
        props = properties(Outputs);

        for i = 1:length(props)
            obj.(props{i}).value = NaN;
            obj.(props{i}).standard_error = NaN;
            obj.(props{i}).unit = "";
        end

    end

    function add(obj, name, value, unit, standard_error)
        % Method to add a new output to the class
        % name:
        % value:
        % standard_error:
        % unit: unit as a string

        if nargin < 5
            standard_error = NaN;
        end

        if isprop(obj, name)
            obj.(name).value = value;
            obj.(name).standard_error = standard_error;
            obj.(name).unit = unit;
        else
            error('Property %s does not exist in Outputs class', name);
        end

    end

    function writeJson(obj, path)
        props = properties(Outputs);
        data = struct();

        for i = 1:length(props)
            data.(props{i}) = obj.(props{i}).value;
        end

        for i = 1:length(props)
            data.(strcat(props{i}, "_ste")) = obj.(props{i}).standard_error;
        end

        for i = 1:length(props)
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

        props = properties(Outputs);

        [dir, name, ~] = fileparts(path);
        path = fullfile(dir, strcat(name, ".h5"));

        if isfile(path) % clear before rewriting
            delete(path)
        end

        for i = 1:length(props)
            h5create(path, strcat("/", props{i}), size(obj.(props{i}).value));
            h5write(path, strcat("/", props{i}), obj.(props{i}).value);
        end

        for i = 1:length(props)
            h5create(path, strcat("/", props{i}, "_ste"), size(obj.(props{i}).standard_error));
            h5write(path, strcat("/", props{i}, "_ste"), obj.(props{i}).standard_error);
        end

        for i = 1:length(props)
            h5create(path, strcat("/", props{i}, "_unit"), [1 1], Datatype = "string");
            h5write(path, strcat("/", props{i}, "_unit"), string(obj.(props{i}).unit));
        end

    end

end

end
