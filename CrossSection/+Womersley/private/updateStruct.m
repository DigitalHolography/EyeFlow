function base = updateStruct(base, update)
    arguments
        base struct
        update struct
    end

    fields = fieldnames(update);
    
    for i = 1:numel(fields)
        f = fields{i};

        if ~isfield(base, f)
            warning_s("Field %s is not present inside base", f);
        end

        base.(f) = update.(f);
    end
end