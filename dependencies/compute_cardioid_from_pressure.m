function [S_card, s_card] = compute_cardioid_from_pressure(sampling_points_inner, sampling_points_outer, pressure_inner, pressure_outer, fs, c)

% prepend a few zeros
%pressure_inner = [zeros(10, size(pressure_inner, 2)); pressure_inner];
%pressure_outer = [zeros(10, size(pressure_outer, 2)); pressure_outer];

f    = linspace(0, fs/2, size(pressure_outer, 1)/2+1).';
f(1) = 1e-20;
k    = 2*pi*f/c;

S_o = fft(pressure_outer);
S_o = S_o(1:end/2+1, :);
S_i = fft(pressure_inner);
S_i = S_i(1:end/2+1, :);

% -------------------------- pressure gradient ----------------------------
dS = (S_o - S_i) ./ repmat(vecnorm(sampling_points_outer - sampling_points_inner), [size(S_o, 1), 1]);
dS = 1./(1i.*k) .* dS;

% % highpass to clean up DC
% win_length = 2;
% win = hann(2*win_length);
% win = win(1:win_length);
% 
% dS(1:win_length, :) = dS(1:win_length, :) .* repmat(win, 1, size(dS, 2));

% fix DC
dS(1, :) = 0;

% 
% % Nyquist bin
% %S_card(end, :) = S_o(end, :) + real(1./(1i.*k(end))) .* real(dS(end, :));  
% 
% % fix DC
% %dS(1, :) = real(dS(2, :));
% %dS(1, :) = 0;
% dS(1, :) = real(dS(1, :));

ds = ifft([dS; conj(flipud(dS(2:end-1, :)))], [], 1, 'symmetric');

% % highpass to clean up DC
% [z, p, k] = butter(9, 100/(fs/2), "high");
% sos = zp2sos(z,p,k);
% %fvtool(sos);
% %h = impz(sos);
% ds = sosfilt(sos, ds);

% % --- fade out ---
% win_length = 500;
% win = hann(2*win_length);
% win = win(win_length+1:end);
% 
% ds(end-win_length+1:end, :) = dS(end-win_length+1:end, :) .* repmat(win, 1, size(dS, 2));

% subtract some DC
%s_card = s_card - repmat(s_card(1, :), size(s_card, 1), 1);

%s_card = s_card(11:end, :);


% cardioid
s_card = pressure_outer + ds;

S_card = 0;

end





% % pressure gradient
% ds = (pressure_outer - pressure_inner) ./ repmat(vecnorm(sampling_points_outer - sampling_points_inner), [size(pressure_outer, 1), 1]);
% 
% % remove DC from ds
% dc = sum(ds, 1);
% 
% ds = ds - repmat(dc, size(ds, 1), 1);
% 
% % integrate ds over time (Riemann sum)
% integral = cumsum(ds, 1) * 1/fs;
% 
% integral = integral - repmat(integral(end, :), size(integral, 1), 1);
% 
% % highpass to clean up DC
% [z, p, k] = butter(9, 100/(fs/2), "high");
% sos = zp2sos(z,p,k);
% %freqz(sos);
% integral = sosfilt(sos, integral);
% 
% % cardioid
% s_card = pressure_outer + c .* integral; 
% 
% S_card = s_card;
% 
% end


