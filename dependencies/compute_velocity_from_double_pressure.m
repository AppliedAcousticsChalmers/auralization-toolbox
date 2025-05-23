function [v] = compute_velocity_from_double_pressure(sampling_points_inner, sampling_points_outer, pressure_inner, pressure_outer, fs, rho)

f     = linspace(0, fs/2, size(pressure_outer, 1)/2+1).';
f(1)  = 1e-20;
omega = 2*pi*f;

S_o = fft(pressure_outer);
S_o = S_o(1:end/2+1, :);
S_i = fft(pressure_inner);
S_i = S_i(1:end/2+1, :);

% pressure gradient
dS = (S_o - S_i) ./ repmat(vecnorm(sampling_points_outer - sampling_points_inner), [size(S_o, 1), 1]);

% velocity
V = -1./(1i.*rho.*omega) .* dS; 

% fix DC
%V(1, :) = abs(V(2, :));
%V(1, :) = real(V(2, :));
V(1, :) = 0;

v = ifft([V; conj(flipud(V(2:end-1, :)))], [], 1, 'symmetric');

end
