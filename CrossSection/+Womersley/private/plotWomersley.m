function plotWomersley(alpha)
% PLOTWOMERSLEY Generates the Womersley velocity profile for a given alpha.
%   plotWomersley(alpha) plots the real part (blue) and the imaginary
%   part (red) of the normalized Womersley curve across the pipe diameter.
%
%   Example: plotWomersley(5)

    % 1. Define radial coordinates (normalized from -1 to 1 for full pipe)
    r_R = linspace(-1, 1, 500);
    
    % 2. Calculate the complex argument for the Bessel function
    % i^(3/2) is equivalent to exp(i * 3*pi/4)
    arg_complex = alpha * exp(1i * 3 * pi / 4);
    
    % 3. Compute the Womersley profile formula
    % u(r) = 1 - J0(alpha * i^(3/2) * r/R) / J0(alpha * i^(3/2))
    J0_wall = besselj(0, arg_complex);
    J0_r = besselj(0, arg_complex * abs(r_R));
    
    u = 1 - (J0_r ./ J0_wall);
    
    % 4. Plotting
    figure('Color', 'w');
    hold on; grid on;
    
    % Real part in blue
    plot(r_R, real(u), 'b', 'LineWidth', 2, 'DisplayName', 'Real Part (In-phase)');
    
    % Imaginary part in red
    plot(r_R, imag(u), 'r', 'LineWidth', 2, 'DisplayName', 'Complex/Imaginary Part (Out-of-phase)');
    
    % Formatting the plot
    xlabel('Normalized Radius (r/R)');
    ylabel('Normalized Velocity');
    title(['Womersley Velocity Profile (\alpha = ', num2str(alpha), ')']);
    legend('Location', 'best');
    
    % Add a line for the pipe walls
    xline([-1, 1], '--k', 'Pipe Wall');
    
    hold off;
end