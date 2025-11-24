classdef DimEnumClass
   enumeration
      D0_scalar
      D1_array
      D2_array
      D3_array
      D4_array
      strings
   end

   methods (Static)
      function leng = get_enum_length()
         keys = enumeration("DimEnumClass");
         leng = numel(keys);
      end

   end
end