function [C_nm, condition_number] = get_c_nm_surface_radial(taps, r_pressure, r_velocity, azi, col, fs, N, c, dynamic_range_dB, precision)
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

    if f(bin) > 10000
        dynamic_range_dB_tmp = 0;
    else
        dynamic_range_dB_tmp = dynamic_range_dB;
    end

    % regularize with hard clipping
    one_over_s = 1 ./ max(s, 10^(-dynamic_range_dB_tmp/20) * max(s)); % regularize
    
    C_nm(bin, :, 1:size(Y_nm_cardioid, 2)) = conj(U) * (one_over_s .* V.');% (V./s.') * U';

    %condition_number(bin) = cond(Y_nm_cardioid(:, :));

end

C_nm(1, :, :) = real(C_nm(2, :, :));

fprintf('\n\n');

end

