classdef Outputs
    %Class to hold the main outputs of the retinal flow analysis pipeline
    
    properties
        HeartBeat
        ArterialMeanVelocity
        ArterialMaximumVelocity
        ArterialMinimumVelocity
        VenousMeanVelocity
        VenousMaximumVelocity
        VenousMinimumVelocity
        ArterialMeanVolumeRate
        ArterialMaximumVolumeRate
        ArterialMinimumVolumeRate
        VenousMeanVolumeRate
        VenousMaximumVolumeRate
        VenousMinimumVolumeRate
        
        
    end
    
    methods
        
        function obj = Outputs()
            % Constructor for the class, fills the properties with default values
            props = properties(OutputsClass);
            for i = 1:length(props)
                obj.(props{i}).value = Nan;
                obj.(props{i}).standard_error = Nan;
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
        
        
        function writeJson(obj,path)
            props = properties(OutputsClass);
            data = struct();
            for i = 1:length(props)
                data.(props{i}) = obj.(props{i}).value;
                data.(strcat(props{i},"_ste")) = obj.(props{i}).standard_error;
                data.(strcat(props{i},"_unit")) = obj.(props{i}).unit;
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