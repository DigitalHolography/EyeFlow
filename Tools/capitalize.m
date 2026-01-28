function res = capitalize(str)
    arguments
        str
    end
    
    res = char(str);
    res(1) = upper(res(1));
    res = string(res);
end