function [C_nm, condition_number] = get_c_nm_surface_cartesian(taps, grid, fs, N, c, max_attenuation, precision)
% Evaluates the gradient in Cartesian coordinates.
%
% All computations are performed in double precision. The result is stored
% in C_nm with precision 'precision'.

%reg_parameter = 1e-3;

f    = linspace(0, fs/2, taps/2+1).';
f(1) = 1e-20; % to avoid kr = 0

k = 2*pi*f/c;

% convert outer sampling points to spherical coordinates
%x = grid.sampling_points_outer(1, :).';
%y = grid.sampling_points_outer(2, :).';
%z = grid.sampling_points_outer(3, :).';

[azi, ele, r_outer] = cart2sph(grid.sampling_points_outer(1, :).', grid.sampling_points_outer(2, :).', grid.sampling_points_outer(3, :).');
col = pi/2 - ele;

% convert inner sampling points to spherical coordinates
%[~, ~, r_inner] = cart2sph(grid.sampling_points_inner(1, :).', grid.sampling_points_inner(2, :).', grid.sampling_points_inner(3, :).');

% compute midpoints between sampling points for gradient computation
% x_mid = (grid.sampling_points_outer(1, :).' + grid.sampling_points_inner(1, :).')/2;
% y_mid = (grid.sampling_points_outer(2, :).' + grid.sampling_points_inner(2, :).')/2;
% z_mid = (grid.sampling_points_outer(3, :).' + grid.sampling_points_inner(3, :).')/2;
% 
% [azi_mid, ele_mid, r_mid] = cart2sph(x_mid, y_mid, z_mid);
% col_mid = pi/2 - ele_mid;

% compute xyz compontents of the vector in direction of the derivative
grid.normal_vec = grid.sampling_points_outer - grid.sampling_points_inner;
grid.normal_vec = grid.normal_vec ./ vecnorm(grid.normal_vec);

% TODO: change azi -> azi_outer etc.
azi_mid = azi;
col_mid = col;
r_mid   = r_outer;

% quadrature matrix (frequency x position x mode) 
C_nm = zeros(length(k), length(r_outer), (N+1)^2, precision);

condition_number = zeros(length(k), 1);

display_progress('Computing the quadrature matrix:');

% loop over frequency
for bin = 1 : length(k)
        
    display_progress(bin/length(k));
    
    % (position x mode)
    Y_nm_pressure = zeros(length(r_outer), (N+1)^2);
    Y_nm_gradient = zeros(length(r_outer), (N+1)^2);
    Y_nm_cardioid = zeros(length(r_outer), (N+1)^2);
    
    % create SH matrices to be inverted
    for n = 0 : N
        
        bessel = sphbesselj(n, k(bin).*r_outer);
         
        for m = -n : n
            
            % pressure
            sh_outer = sphharm(n, m, col, azi, 'real');  
            Y_nm_pressure(:, n^2+n+m+1) = 4*pi .* 1i^n .* bessel .* sh_outer;
            
            % gradient            
            Y_nm_gradient(:, n^2+n+m+1) = Y_nm_gradient(:, n^2+n+m+1) + grid.normal_vec(1, :).' .* 4*pi .* 1i^n .* d_sh_mode(n, m, r_mid, col_mid, azi_mid, k(bin), 'dx');
            Y_nm_gradient(:, n^2+n+m+1) = Y_nm_gradient(:, n^2+n+m+1) + grid.normal_vec(2, :).' .* 4*pi .* 1i^n .* d_sh_mode(n, m, r_mid, col_mid, azi_mid, k(bin), 'dy');
            Y_nm_gradient(:, n^2+n+m+1) = Y_nm_gradient(:, n^2+n+m+1) + grid.normal_vec(3, :).' .* 4*pi .* 1i^n .* d_sh_mode(n, m, r_mid, col_mid, azi_mid, k(bin), 'dz');
            
            % cardioid
            Y_nm_cardioid(:, n^2+n+m+1) = Y_nm_pressure(:, n^2+n+m+1) + 1./(1i.*k(bin)) .* Y_nm_gradient(:, n^2+n+m+1);
                    
            % ------------ regularization by soft clipping ----------------
            
            if n <= 1
                % Relax the regularization for orders 0 and 1 for low
                % frequencies
                alpha = 10^(     80        /20) * ones(size(r_outer));
            else 
                alpha = 10^(max_attenuation/20) * ones(size(r_outer));
            end
            
            if n > 1
                % regularize high orders and low frequencies harder
                alpha(k(bin).*r_outer < n) = 10^(40/20); 
            end
        
            Y_nm_cardioid(  :, n^2+n+m+1) = soft_clipping(Y_nm_cardioid(  :, n^2+n+m+1), alpha);
            
        end
        
    end
    
    % --------------- regularization by soft clipping ---------------------

    % finally, invert the matrix; if C_nm is single precision, the data
    % will be converted automatically from double
    C_nm(bin, :, :) = pinv(Y_nm_cardioid).';
    condition_number(bin) = cond(Y_nm_cardioid);
        
end

fprintf('\n\n');

end

% -------------------------------------------------------------------------
function [data] = soft_clipping(data, alpha)

% avoid dividing by 0
data(data == 0) = 1e-30;    
    
% soft clipping
data_norm = abs(data) ./ data; 
data = 1./((2.*alpha/pi) .* data_norm .* atan(pi./(2 .* alpha .* abs(data))));

end
