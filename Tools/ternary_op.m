function res = ternary_op(cond, T_func, F_func)
    arguments
        cond
        T_func
        F_func
    end

    if cond
        res = T_func();
    else
        res = F_func();
    end
end