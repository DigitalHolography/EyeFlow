function K = calculateKFactor(lambda)
    % Calculates the flow moment factor K(alpha).
    % K(alpha) = 1 - (2*J1(lambda) / (lambda * J0(lambda)))
    
    if lambda == 0 || besselj(0, lambda) == 0
        K = complex(NaN, NaN); % Avoid division by zero
        return;
    end
    
    K = 1 - (2 * besselj(1, lambda)) / (lambda * besselj(0, lambda));
end