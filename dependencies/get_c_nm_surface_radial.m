function [C_nm, condition_number] = get_c_nm_surface_radial(taps, grid_data, fs, N, c, dynamic_range_dB, precision)
% Evaluates the gradient in radial direction.
%
% All computations are performed in double precision. The result is stored
% in C_nm with precision 'precision'.

% maximum magnitude that occurs in the matrix to be inverted; determines
% the regularization limit
max_magnitude = -inf;

f    = linspace(0, fs/2, taps/2+1).';
f(1) = 1e-20; % to avoid kr = 0

k = 2*pi*f/c;

% for measuring the pressure
[azi, ele, r_pressure] = cart2sph(grid_data.sampling_points_outer(1, :).', grid_data.sampling_points_outer(2, :).', grid_data.sampling_points_outer(3, :).');
col = pi/2 - ele;

[~, ~, r_gradient] = cart2sph(grid_data.sampling_points_inner(1, :).', grid_data.sampling_points_inner(2, :).', grid_data.sampling_points_inner(3, :).'); 
r_gradient = (r_gradient + r_pressure)/2;

% quadrature matrix (frequency x position x mode) 
C_nm = zeros(length(k), length(r_pressure), (N+1)^2, precision);

condition_number = zeros(length(k), 1);

display_progress('Computing the quadrature matrix:');

% (position x mode x bin)
Y_nm_cardioid = zeros(length(r_pressure), (N+1)^2, length(k));
    
% ----------------- assemble the matrix to be inverted --------------------
for bin = 1 : length(k)
        
    display_progress(bin/length(k));
    
    % create SH matrices to be inverted
    for n = 0 : N
        
        bessel       = sphbesselj(n, k(bin).*r_pressure);
        bessel_prime = 1/(2*n+1) * (n * sphbesselj(n-1, k(bin).*r_gradient) - (n+1) * sphbesselj(n+1, k(bin).*r_gradient));
         
        for m = -n : n
            
            % cardioid
            Y_nm_cardioid(:, n^2+n+m+1, bin) = 4*pi .* 1i^n .* (bessel - 1i .* bessel_prime) .* sphharm(n, m, col, azi, 'real');
             
        end
        
    end
    
end

% ------------------- regularize and invert the matrix --------------------
for bin = 1 : length(k)
    
    Y_nm_cardioid(:, :, bin) = soft_clipping(Y_nm_cardioid(:, :, bin), dynamic_range_dB);
    
    % finally, invert the matrix; if C_nm is single precision, the data
    % will be converted automatically from double
    C_nm(bin, :, :) = pinv(Y_nm_cardioid(:, :, bin)).';
    
    condition_number(bin) = cond(Y_nm_cardioid(:, :, bin));
        
end

fprintf('\n\n');

end

% -------------------------------------------------------------------------
function [data] = soft_clipping(data, dynamic_range_dB)

alpha = 10^(dynamic_range_dB/20);

% avoid dividing by 0
data(data == 0) = 1e-30;
data(abs(data) < 1e-30) = data(abs(data) < 1e-30) ./ abs(data(abs(data) < 1e-30)) * 1e-30; 

maximum = max(abs(data(:)));
minimum = min(abs(data(:)));

%disp(20*log10(abs([minimum(1), maximum(1)])));

if maximum(1)/minimum(1) > alpha
    % alpha is the reciprocal value of the smallest permitted magnitude
    alpha = alpha/maximum(1);
end

alpha = alpha * ones(size(data));

% soft clipping
data_norm = abs(data) ./ data; 
data = 1./((2.*alpha/pi) .* data_norm .* atan(pi./(2 .* alpha .* abs(data))));

end






