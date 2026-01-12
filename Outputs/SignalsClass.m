classdef SignalsClass < handle
% Class to hold the main output signals of the retinal flow analysis pipeline

properties
    ArterialVelocity
    VenousVelocity
    ArterialVolumeRate
    VenousVolumeRate
    % ArterialVelocityInterpolated
    % VenousVelocityInterpolated
    TransFunctionModLog10
    TransFunctionPhaseDegrees
    % ArterialVolumeRateInterpolated
    % VenousVolumeRateInterpolated

end

methods

    function obj = SignalsClass()
        % Constructor for the class, fills the properties with default values
        props = properties(obj);

        for i = 1:length(props)
            obj.(props{i}).yvalues = NaN;
            obj.(props{i}).xvalues = NaN;
            obj.(props{i}).ystandard_errors = NaN;
            obj.(props{i}).xunit = "";
            obj.(props{i}).yunit = "";
        end

    end

    function add(obj, name, yvalues, yunit, xvalues, xunit, ystandard_errors)
        % Method to add a new output to the class
        % name:
        % yvalues:
        % yunit: unit as a string
        % xvalues:
        % xunit: unit as a string
        % ystandard_errors:

        if nargin < 7
            ystandard_errors = NaN;
        end

        if isprop(obj, name)
            obj.(name).yvalues = yvalues;
            obj.(name).xvalues = xvalues;
            obj.(name).ystandard_errors = ystandard_errors;
            obj.(name).xunit = xunit;
            obj.(name).yunit = yunit;
        else
            error('Property %s does not exist in Signals class', name);
        end

    end

    function writeHdf5(obj, path)

        props = properties(obj);

        for i = 1:length(props)
            h5create(path, strcat("/", "Signals", "/", props{i}, "/", props{i}, "_y"), size(obj.(props{i}).yvalues));
            h5write(path, strcat("/", "Signals", "/", props{i}, "/", props{i}, "_y"), obj.(props{i}).yvalues);
            h5writeatt(path, strcat("/", "Signals", "/", props{i}, "/", props{i}, "_y"), "unit", obj.(props{i}).yunit);
        end

        for i = 1:length(props)
            h5create(path, strcat("/", "Signals", "/", props{i}, "/", props{i}, "_y_ste"), size(obj.(props{i}).ystandard_errors));
            h5write(path, strcat("/", "Signals", "/", props{i}, "/", props{i}, "_y_ste"), obj.(props{i}).ystandard_errors);
            h5writeatt(path, strcat("/", "Signals", "/", props{i}, "/", props{i}, "_y_ste"), "unit", obj.(props{i}).yunit);
        end

        for i = 1:length(props)
            h5create(path, strcat("/", "Signals", "/", props{i}, "/", props{i}, "_x"), size(obj.(props{i}).xvalues));
            h5write(path, strcat("/", "Signals", "/", props{i}, "/", props{i}, "_x"), obj.(props{i}).xvalues);
            h5writeatt(path, strcat("/", "Signals", "/", props{i}, "/", props{i}, "_x"), "unit", obj.(props{i}).xunit);
        end

    end

end

end
