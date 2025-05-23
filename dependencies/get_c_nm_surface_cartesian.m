function [C_nm, condition_number] = get_c_nm_surface_cartesian(taps, r_pressure, azi_pressure, col_pressure, r_velocity, azi_velocity, col_velocity, normal_vector, fs, N, c, dynamic_range_permitted_dB, precision)
% Evaluates the gradient in Cartesian coordinates.
%
% All computations are performed in double precision. The result is stored
% in C_nm with precision 'precision'.

%r_max = max(r_pressure(:));
%r_max = r_max(1);

% maximum magnitude that occurs in the matrix to be inverted; determines
% the regularization limit
%max_magnitude = 1;%4*pi * sphbesselj(0, 0) * sphharm(0, 0, 0, 0, 'real');

f    = linspace(0, fs/2, taps/2+1).';
f(1) = 1e-20; % to avoid kr = 0

k = 2*pi*f/c;

% assure unit length
normal_vector = normal_vector ./ vecnorm(normal_vector);

% quadrature matrix (frequency x position x mode) 
C_nm = zeros(length(k), length(r_pressure), (N+1)^2, precision);

condition_number = zeros(length(k), 1);

display_progress('Computing the quadrature matrix:');

% loop over frequency
for bin = 2 : length(k) % skip DC
        
    display_progress(bin/length(k));
    
    % (position x mode)
    Y_nm_pressure = zeros(length(r_pressure), (N+1)^2);
    Y_nm_gradient = zeros(length(r_pressure), (N+1)^2);
    Y_nm_cardioid = zeros(length(r_pressure), (N+1)^2);
    
    % create SH matrices to be inverted
    for n = 0 : N
        
        bessel = sphbesselj(n, k(bin).*r_pressure);
         
        for m = -n : n
            
            % pressure
            sh_outer = sphharm(n, m, col_pressure, azi_pressure, 'real');  
            Y_nm_pressure(:, n^2+n+m+1) = 4*pi .* 1i^n .* bessel .* sh_outer;
            
            % gradient            
            Y_nm_gradient(:, n^2+n+m+1) = Y_nm_gradient(:, n^2+n+m+1) + normal_vector(1, :).' .* 4*pi .* 1i^n .* d_sh_mode(n, m, r_velocity, col_velocity, azi_velocity, k(bin), 'dx');
            Y_nm_gradient(:, n^2+n+m+1) = Y_nm_gradient(:, n^2+n+m+1) + normal_vector(2, :).' .* 4*pi .* 1i^n .* d_sh_mode(n, m, r_velocity, col_velocity, azi_velocity, k(bin), 'dy');
            Y_nm_gradient(:, n^2+n+m+1) = Y_nm_gradient(:, n^2+n+m+1) + normal_vector(3, :).' .* 4*pi .* 1i^n .* d_sh_mode(n, m, r_velocity, col_velocity, azi_velocity, k(bin), 'dz');
            
            % cardioid
            Y_nm_cardioid(:, n^2+n+m+1) = Y_nm_pressure(:, n^2+n+m+1) + 1./(1i.*k(bin)) .* Y_nm_gradient(:, n^2+n+m+1);
                    
            % ------------ regularization by soft clipping ----------------
            
            %max_magnitude = max(abs(Y_nm_cardioid(:, n^2+n+m+1)));

            %alpha = 10^(dynamic_range_dB/20);
            %alpha = 10^(dynamic_range_dB/20)/max_magnitude * ones(size(r_pressure));

            % if n > 1
            %     % regularize high orders and low frequencies harder
            %     alpha(k(bin).*r_max < n) = 10^(0/20)/max_magnitude;
            % end

            %%%%%%%%%%%Y_nm_cardioid(:, n^2+n+m+1) = soft_clipping(Y_nm_cardioid( :, n^2+n+m+1), alpha);
            
        end
        
    end
      
    %[U, s, V] = svd(Y_nm_cardioid, 'econ', 'vector');
    [U, s, V] = svd(Y_nm_cardioid, 'econ'); s = diag(s); % for older MATLAB versions

    dynamic_range_actual_dB   = 20*log10(s(1)/s(end));
    dynamic_range_to_apply_dB = get_dynamic_range_to_apply_dB(f(bin), dynamic_range_actual_dB, dynamic_range_permitted_dB);

    % reduce rank by 1 for safety
    s = s(1:end-1);
    U = U(:, 1:end-1);
    V = V(:, 1:end-1);

    % invert the matrix
    s = 1 ./ max(s, 10^(-dynamic_range_to_apply_dB/20) * max(s));
    C_nm(bin, :, :) = conj(U) * (s .* V.');% (V./s.') * U';

    %condition_number(bin) = 1; %cond(Y_nm_cardioid);

end

% fix DC
C_nm(1, :, :) = real(C_nm(2, :, :)); % this avoids a high condition number at DC

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

else
    
    if dynamic_range_actual_dB < dynamic_range_permitted_dB(1) + head_room_dB
        dynamic_range_to_apply_dB = dynamic_range_actual_dB - head_room_dB;
    else
        dynamic_range_to_apply_dB = dynamic_range_permitted_dB(1);
    end

end

dynamic_range_to_apply_dB = max(0, dynamic_range_to_apply_dB);

end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% function [data] = soft_clipping(data, alpha)
% 
% % avoid dividing by 0
% data(data == 0) = 1e-30;    
% 
% % soft clipping
% data_norm = abs(data) ./ data; 
% data = 1./((2.*alpha/pi) .* data_norm .* atan(pi./(2 .* alpha .* abs(data))));
% 
% end
