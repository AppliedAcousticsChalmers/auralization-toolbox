function [C_nm, condition_number] = get_c_nm_surface_radial(taps, r_pressure, r_velocity, azi, col, fs, N, c, dynamic_range_permitted_dB, precision)
% Evaluates the gradient in radial direction.
% 
% r_pressure: radius at which the pressure is observed 
%
% r_velocity: radius at which the velocity is observed (for double layers,
%             this one is slightly different than r_pressure).
%
% All computations are performed in double precision. The result is stored
% in C_nm with precision 'precision'.

f    = linspace(0, fs/2, taps/2+1).';
f(1) = 1e-20; % to avoid kr = 0

k = 2*pi*f/c;

% quadrature matrix (frequency x position x mode) 
C_nm = zeros(length(k), length(azi), (N+1)^2, precision);

condition_number = ones(length(k), 1);

display_progress('Computing the quadrature matrix:');
    
% ----------------- assemble the matrix to be inverted --------------------
for bin = 2 : length(k) % DC is badly conditioned
        
    display_progress(bin/length(k));
    
    % create SH matrix to be inverted (position x mode)
    Y_nm_cardioid = zeros(length(azi), (N+1)^2);
    
    for n = 0 : N
        
        bessel       = sphbesselj(n, k(bin).*r_pressure);
        bessel_prime = 1/(2*n+1) * (n * sphbesselj(n-1, k(bin).*r_velocity) - (n+1) * sphbesselj(n+1, k(bin).*r_velocity));
         
        for m = -n : n
            
            % cardioid
            Y_nm_cardioid(:, n^2+n+m+1) = 4*pi .* 1i^n .* (bessel - 1i .* bessel_prime) .* sphharm(n, m, col, azi, 'real');
              
        end
        
    end

    % ----------------- regularize and invert the matrix ------------------

    %[U, s, V] = svd(Y_nm_cardioid, 'econ', 'vector');
    [U, s, V] = svd(Y_nm_cardioid, 'econ'); s = diag(s); % for older MATLAB versions

    % invert
    %one_over_s = 1./s;

    % regularize with soft clipping
    %one_over_s = soft_clip_sv(one_over_s, dynamic_range_dB);

    dynamic_range_actual_dB   = 20*log10(s(1)/s(end));
    dynamic_range_to_apply_dB = get_dynamic_range_to_apply_dB(f(bin), dynamic_range_actual_dB, dynamic_range_permitted_dB);

    % reduce rank by 1 for safety
    s = s(1:end-1);
    U = U(:, 1:end-1);
    V = V(:, 1:end-1);

    % invert the matrix
    one_over_s = 1 ./ max(s, 10^(-dynamic_range_to_apply_dB/20) * max(s)); % regularize
    
    C_nm(bin, :, 1:size(Y_nm_cardioid, 2)) = conj(U) * (one_over_s .* V.');% (V./s.') * U';

    %condition_number(bin) = cond(Y_nm_cardioid(:, :));

end

C_nm(1, :, :) = real(C_nm(2, :, :));

fprintf('\n\n');

end

% -------------------------------------------------------------------------
function dynamic_range_to_apply_dB = get_dynamic_range_to_apply_dB(f, dynamic_range_actual_dB, dynamic_range_permitted_dB)

head_room_dB = 3; % make sure that some amount of regularization is always being applied

% --- determine how much to regularize ---
if f > 10000 % 10 kHz

    if dynamic_range_actual_dB < dynamic_range_permitted_dB(3) + head_room_dB
        dynamic_range_to_apply_dB = dynamic_range_actual_dB - head_room_dB;
    else
        dynamic_range_to_apply_dB = dynamic_range_permitted_dB(3);
    end

elseif f > 200 % 200 Hz - 10 kHz

    if dynamic_range_actual_dB < dynamic_range_permitted_dB(2) + head_room_dB
        dynamic_range_to_apply_dB = dynamic_range_actual_dB - head_room_dB;
    else
        dynamic_range_to_apply_dB = dynamic_range_permitted_dB(2);
    end

else % f < 200 Hz
    
    if dynamic_range_actual_dB < dynamic_range_permitted_dB(1) + head_room_dB
        dynamic_range_to_apply_dB = dynamic_range_actual_dB - head_room_dB;
    else
        dynamic_range_to_apply_dB = dynamic_range_permitted_dB(1);
    end

end

dynamic_range_to_apply_dB = max(0, dynamic_range_to_apply_dB);

end
% -------------------------------------------------------------------------
