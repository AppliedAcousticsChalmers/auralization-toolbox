function [S_card, s_card] = compute_cardioid_from_pressure(sampling_points_inner, sampling_points_outer, pressure_inner, pressure_outer, fs, c)

f    = linspace(0, fs/2, size(pressure_outer, 1)/2+1).';
f(1) = 1e-20;
k    = 2*pi*f/c;

S_o = fft(pressure_outer);
S_o = S_o(1:end/2+1, :);
S_i = fft(pressure_inner);
S_i = S_i(1:end/2+1, :);

% pressure gradient
dS = (S_o - S_i) ./ repmat(vecnorm(sampling_points_outer - sampling_points_inner), [size(S_o, 1), 1]);

% cardioid
S_card = S_o + 1./(1i.*k) .* dS; 

% fix DC
S_card(1, :) = abs(S_card(2, :));

s_card = ifft([S_card; conj(flipud(S_card(2:end-1, :)))], [], 1, 'symmetric');

end