classdef DimNameEnumClass
    % Add all possible Attributes values inside
    enumeration
        % Complex values
        RealPart
        ImagPart

        BranchIdx
        CircleIdx
        SegmentIdx

        Harmonic
    end

    methods (Static)
        function leng = get_enum_length()
            keys = enumeration("DimNameEnumClass");
            leng = numel(keys);
        end
    end
end