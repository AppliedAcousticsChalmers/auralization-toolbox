function [C_nm, condition_number] = get_c_nm_volumetric(taps, grid_data, fs, N, c, max_attenuation, sphharm_type, precision)
% All computations are performed in double precision. The result is stored
% in C_nm with precision 'precision'.
% 
% C_nm: Quadrature matrix in frequency domain

f    = linspace(0, fs/2, taps/2+1).';
f(1) = 1e-20; % to avoid kr = 0

k = 2*pi*f/c;

% convert grid to spherical coordinates
[azi, ele, r] = cart2sph(grid_data.sampling_points(1, :).', grid_data.sampling_points(2, :).', grid_data.sampling_points(3, :).');
col = pi/2 - ele;

% quadrature matrix
C_nm             = zeros(length(k), length(r), (N+1)^2, precision);
condition_number = zeros(length(k), 1);

display_progress('Computing the quadrature matrix:');

% loop over frequency
for bin = 1 : length(k)
        
    display_progress(bin/length(k));
    
    Y_nm = zeros(length(r), (N+1)^2);
                
    % create SH matrix to be inverted
    for n = 0 : N
        
        bessel_term = 4*pi .* 1i^n .* sphbesselj(n, k(bin).*r);
        
        for m = -n : n
            
            %Y_nm_tmp = bessel_term .* sphharm(n, m, col, azi, sphharm_type);            
            %Y_nm(:, n^2+n+m+1) = Y_nm_tmp;
            
            Y_nm_tmp = bessel_term .* sphharm(n, m, col, azi, sphharm_type);
            
            alpha = 10^(max_attenuation/20) * ones(size(r));
            % regularize high orders and low frequencies harder
            alpha(k(bin).*r < n) = 10^(30/20); 
            
            % avoid dividing by 0
            Y_nm_tmp(Y_nm_tmp == 0) = 1e-30;

            Y_nm_norm = abs(Y_nm_tmp) ./ Y_nm_tmp; 
            Y_nm_tmp = 1./((2.*alpha/pi) .* Y_nm_norm .* atan(pi./(2 .* alpha .* abs(Y_nm_tmp))));
            
            Y_nm(:, n^2+n+m+1) = Y_nm_tmp;
        end
        
    end
        
%     % set the regularization parameter
%     %max_Y_nm = max(abs(Y_nm(:)));
%     alpha = 10^(dynamic_range_dB/20) * ones(size(Y_nm));% / max_Y_nm(1);
% 
%     % ----- optional ----
%     for n = 1 : N
%         for m = -n : n
%             % regularize high orders and low frequencies harder
%             alpha(k(bin).*r < n, n^2+n+m+1) = 10^(40/20);% / max_Y_nm(1);
%         end
%     end
% 
%     % avoid dividing by 0
%     Y_nm(abs(Y_nm) < 1e-30) = 1e-30;
% 
%     % regularization by soft clipping
%     Y_nm_norm = abs(Y_nm) ./ Y_nm; 
%     Y_nm_tmp = 1./((2.*alpha/pi) .* Y_nm_norm .* atan(pi./(2 .* alpha .* abs(Y_nm_tmp))));

    % finally, invert the matrix; if C_nm is single precision, the data
    % will be converted automatically from double
    C_nm(bin, :, :) = pinv(Y_nm).';
    
    condition_number(bin) = cond(Y_nm);
        
end

fprintf('\n\n');


