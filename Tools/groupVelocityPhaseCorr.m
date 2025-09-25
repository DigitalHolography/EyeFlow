function vg = groupVelocityPhaseCorr(U, x, t)
% groupVelocityPhaseCorr Estimate group velocity from U(x,t) using phase correlation
%
%   vg = groupVelocityPhaseCorr(U, x, t)
%
%   INPUTS:
%     U : 2D matrix (nx × nt), field U(x,t)
%         rows = spatial points (x), columns = time snapshots (t)
%     x : vector of spatial coordinates (nx × 1)
%     t : vector of temporal coordinates (1 × nt)
%
%   OUTPUT:
%     vg : estimated group velocity (same units as x/t)

% --- Check inputs ---
[nx, nt] = size(U);

% --- FFT along space and time ---
Uxw = fftshift(fftn(U));    % 2D FFT: kx × omega
Sk = abs(Uxw).^2;           % power spectrum

% --- Frequency axes ---
dx = mean(diff(x));
dt = mean(diff(t));
kx = fftshift(fftspace(nx, dx));     % spatial frequency (rad/unit)
w  = fftshift(fftspace(nt, dt));     % temporal frequency (rad/unit time)

% --- Find spectral peak ---
[~, idx] = max(Sk(:));
[ik, iw] = ind2sub(size(Sk), idx);

% --- Extract k and omega at peak ---
k_peak = kx(ik);
w_peak = w(iw);

% --- Group velocity estimate ---
vg = w_peak / k_peak;

end

% --- Helper: FFT frequency axis in radians ---
function f = fftspace(n, d)
    % cycles per unit
    freq = (-floor(n/2):ceil(n/2)-1)/(n*d);
    % radians per unit
    f = 2*pi*freq;
end