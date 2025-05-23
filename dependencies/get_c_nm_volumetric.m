function [C_nm, condition_number] = get_c_nm_volumetric(taps, sampling_points, fs, N, c, dynamic_range_permitted_dB, sphharm_type, precision)
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

    dynamic_range_actual_dB   = 20*log10(s(1)/s(end));
    dynamic_range_to_apply_dB = get_dynamic_range_to_apply_dB(f(bin), dynamic_range_actual_dB, dynamic_range_permitted_dB);

    % reduce rank by 1 for safety
    s = s(1:end-1);
    U = U(:, 1:end-1);
    V = V(:, 1:end-1);

    % invert the matrix
    s = 1 ./ max(s, 10^(-dynamic_range_to_apply_dB/20) * max(s)); % regularize
    C_nm(bin, :, :) = conj(U) * (s .* V.');% (V./s.') * U';
    
    %condition_number(bin) = cond(Y_nm);
        
end

% fix DC
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

