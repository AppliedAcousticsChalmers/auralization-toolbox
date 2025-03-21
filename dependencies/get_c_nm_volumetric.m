function [C_nm, condition_number] = get_c_nm_volumetric(taps, sampling_points, fs, N, c, dynamic_range_dB, sphharm_type, precision)
% All computations are performed in double precision. The result is stored
% in C_nm with precision 'precision'.
% 
% C_nm: Quadrature matrix in frequency domain

f    = linspace(0, fs/2, taps/2+1).';
f(1) = 1e-20; % to avoid kr = 0

k = 2*pi*f/c;

% convert grid to spherical coordinates
[azi, ele, r] = cart2sph(sampling_points(1, :).', sampling_points(2, :).', sampling_points(3, :).');
col = pi/2 - ele;

% quadrature matrix
C_nm             = zeros(length(k), length(r), (N+1)^2, precision);
condition_number = zeros(length(k), 1);

display_progress('Computing the quadrature matrix:');

% for regularizing the dynamic range
%max_magnitude = 4*pi .* 1i^0 .* sphbesselj(0, 0) .* sphharm(0, 0, 0, 0, sphharm_type);

% loop over frequency (DC can be very badly conditioned)
for bin = 2 : length(k)
        
    display_progress(bin/length(k));
    
    Y_nm = zeros(length(r), (N+1)^2);
                
    % create SH matrix to be inverted
    for n = 0 : N
        
        bessel_term = 4*pi .* 1i^n .* sphbesselj(n, k(bin).*r);

        %max_magnitude = max(abs(bessel_term(:)));
       
        for m = -n : n
            
            %Y_nm_tmp = bessel_term .* sphharm(n, m, col, azi, sphharm_type);            
            %Y_nm(:, n^2+n+m+1) = Y_nm_tmp;
            
            Y_nm_tmp = bessel_term .* sphharm(n, m, col, azi, sphharm_type);

            % alpha = 10^(dynamic_range_dB/20)/max_magnitude * ones(size(r));
            % % regularize high orders and low frequencies harder
            % %alpha(k(bin).*r < n) = 10^(40/20)/max_magnitude; 
            % 
            % % avoid dividing by 0
            % Y_nm_tmp(Y_nm_tmp == 0) = 1e-30;
            % 
            % Y_nm_norm = abs(Y_nm_tmp) ./ Y_nm_tmp; 
            % Y_nm_tmp = 1./((2.*alpha/pi) .* Y_nm_norm .* atan(pi./(2 .* alpha .* abs(Y_nm_tmp))));

            Y_nm(:, n^2+n+m+1) = Y_nm_tmp;
        end
        
    end
    
    % finally, invert the matrix; if C_nm is single precision, the data
    % will be converted automatically from double
    %C_nm(bin, :, :) = pinv(Y_nm).';

    %[U, s, V] = svd(Y_nm_cardioid, 'econ', 'vector');
    [U, s, V] = svd(Y_nm, 'econ'); s = diag(s); % for older MATLAB versions

    if f(bin) > 10000 % 10 kHz
        dynamic_range_dB_tmp = 0;
    else
        dynamic_range_dB_tmp = dynamic_range_dB;
    end

    s = 1 ./ max(s, 10^(-dynamic_range_dB_tmp/20) * max(s)); % regularize
    C_nm(bin, :, :) = conj(U) * (s .* V.');% (V./s.') * U';
    
    %condition_number(bin) = cond(Y_nm);
        
end

% fix DC
C_nm(1, :, :) = real(C_nm(2, :, :));

fprintf('\n\n');


