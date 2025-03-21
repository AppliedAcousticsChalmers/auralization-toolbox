function [C_nm, condition_number] = get_c_nm_surface_cartesian(taps, r_pressure, azi_pressure, col_pressure, r_velocity, azi_velocity, col_velocity, normal_vector, fs, N, c, dynamic_range_dB, precision)
% Evaluates the gradient in Cartesian coordinates.
%
% All computations are performed in double precision. The result is stored
% in C_nm with precision 'precision'.

r_max = max(r_pressure(:));
r_max = r_max(1);

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
    
    % --------------- regularization by soft clipping ---------------------

    % max_magnitude = max(abs(Y_nm_cardioid(:)));
    % 
    % %if n <= 1
    %     % Relax the regularization for orders 0 and 1 for low
    %     % frequencies
    %     %alpha = 10^(60/20);
    % %    alpha = 10^(      60        /20)/max_magnitude * ones(size(r_pressure));
    % %else 
    %     %alpha = 10^(dynamic_range_dB/20);
    % dynamic_range_dB = 60;
    %     alpha = 10^(dynamic_range_dB/20)/max_magnitude * ones(size(Y_nm_cardioid));
    % % %end
    % 
    % if n > 1
    %    % regularize high orders and low frequencies harder
    %    alpha(k(bin).*r_pressure < n) = 10^(20/20); 
    % end
    % 
    % Y_nm_cardioid = soft_clipping(Y_nm_cardioid, alpha);

    % finally, invert the matrix; if C_nm is single precision, the data
    % will be converted automatically from double
    %%%%%%%%%%%%% C_nm(bin, :, :) = pinv(Y_nm_cardioid, 1e-1).';
    
    %[U, s, V] = svd(Y_nm_cardioid, 'econ', 'vector');
    [U, s, V] = svd(Y_nm_cardioid, 'econ'); s = diag(s); % for older MATLAB versions

    if f(bin) > 10000 % 10 kHz
        dynamic_range_dB_tmp = 0;
    else
        dynamic_range_dB_tmp = dynamic_range_dB;
    end

    s = 1 ./ max(s, 10^(-dynamic_range_dB_tmp/20) * max(s)); % regularize
    C_nm(bin, :, :) = conj(U) * (s .* V.');% (V./s.') * U';

    condition_number(bin) = 1; %cond(Y_nm_cardioid);

end

% fix DC
C_nm(1, :, :) = real(C_nm(2, :, :)); % this avoids a high condition number at DC

fprintf('\n\n');

end

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
