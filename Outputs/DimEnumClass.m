classdef DimEnumClass
% All possible "Structures" that can be saved inside h5 (in their own folder)
% To add one, you must update:
%     - The enumeration
%     - The validate function

properties (SetAccess = immutable)
    rank double
end

enumeration
    Array0D_FrequencyAnalysis (0)
    Array1D_FrequencyAnalysis (1)
    Array2D_FrequencyAnalysis (2)
    Array3D_FrequencyAnalysis (3)
    Array4D_FrequencyAnalysis (4)
    Array5D_FrequencyAnalysis (5)
    Strings (1)
    Other (NaN)
end

methods

    function obj = DimEnumClass(rank)

        arguments
            rank double
        end

        obj.rank = rank;
    end

    function res_bool = validate(obj, data)

        arguments
            obj
            data
        end

        if obj == DimEnumClass.Other
            res_bool = true;
            return;
        end

        % strings
        if isstring(data) || ischar(data) || iscellstr(data)
            res_bool = obj == DimEnumClass.Strings;
            return;
        end

        % ndims gives 2 for a scalar and a vector
        if obj.rank == 0
            res_bool = isscalar(data);
        elseif obj.rank == 1
            res_bool = isvector(data) && ~isscalar(data);
        else
            res_bool = ndims(data) == obj.rank;
        end

    end

end

methods (Static)

    function leng = get_enum_length()
        keys = enumeration("DimEnumClass");
        leng = numel(keys);
    end

    function obj = enum_from(data)
        keys = enumeration("DimEnumClass");
        leng = numel(keys);

        for i = 1:leng
            cur_enum = keys(i);

            if cur_enum.validate(data)
                obj = cur_enum;
                return
            end

        end

        warning("DimEnumClass:enum_from", "The Data inserted has found no valid match.")
    end

end

end
