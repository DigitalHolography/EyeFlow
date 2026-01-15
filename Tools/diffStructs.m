function [isSame, diffMsg] = diffStructs(s1, s2)
    isSame = true;
    diffMsg = "";

    if isstruct(s1) ~= isstruct(s2)
        isSame = false;
        diffMsg = sprintf('Type mismatch: One is struct, one is value.');
        return;
    end

    f1 = sort(fieldnames(s1));
    f2 = sort(fieldnames(s2));

    if ~isequal(f1, f2)
        isSame = false;

        missingInS2 = setdiff(f1, f2);
        missingInS1 = setdiff(f2, f1);

        diffMsg = sprintf("Mismatch fields. \nIn S1 only: %s \nIn S2 only: %s", ...
            strjoin(missingInS2, ', '), strjoin(missingInS1, ', '));
        return;
    end

    % Nested strcuts
    for i = 1:numel(f1)
        field = f1{i};
        
        val1 = s1(1).(field);
        val2 = s2(1).(field);
        
        if isstruct(val1) && isstruct(val2)
            [subSame, subMsg] = diffStructs(val1, val2);
            if ~subSame
                isSame = false;
                diffMsg = sprintf('Mismatch in sub-struct ".%s": %s', field, subMsg);
                return;
            end
        elseif isstruct(val1) ~= isstruct(val2)
            isSame = false;
            diffMsg = sprintf('Type mismatch in field ".%s": One is struct, one is value.', field);
            return;
        end
    end
end