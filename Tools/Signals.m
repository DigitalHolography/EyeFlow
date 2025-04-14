classdef Signals < handle
    %Class to hold the main output signals of the retinal flow analysis pipeline
    
    properties
        ArterialVelocity
        VenousVelocity
        % ArterialVelocityInterpolated
        % VenousVelocityInterpolated
        ArterialVolumeRate
        VenousVolumeRate
        % ArterialVolumeRateInterpolated
        % VenousVolumeRateInterpolated
    end
    
    methods
        
        function obj = Signals()
        end
        
        function initSignals(obj)
            % Constructor for the class, fills the properties with default values
            props = properties(Signals);
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
        
        
        function writeJson(obj,path)
            props = properties(Signals);
            data = struct();
            for i = 1:length(props)
                data.(strcat(props{i},"_y")) = obj.(props{i}).yvalues;
            end
            for i = 1:length(props)
                data.(strcat(props{i},"_ste")) = obj.(props{i}).ystandard_errors;
            end
            for i = 1:length(props)
                data.(strcat(props{i},"_x")) = obj.(props{i}).xvalues;
            end
            for i = 1:length(props)
                data.(strcat(props{i},"_yunit")) = obj.(props{i}).yunit;
            end
            for i = 1:length(props)
                data.(strcat(props{i},"_xunit")) = obj.(props{i}).xunit;
            end
            jsonText = jsonencode(data,"PrettyPrint",true);
            fid = fopen(path, 'w');
            if fid == -1
                error('Cannot open file for writing: %s', path);
            end
            fwrite(fid, jsonText, 'char');
            fclose(fid);
        end
    end
end