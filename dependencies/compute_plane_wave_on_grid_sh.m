function [pressure, velocity] = compute_plane_wave_on_grid_sh(sampling_points, pw_incidence_azi_deg, pw_incidence_ele_deg, c, taps_pw, fs, grid_shape, normal_vector, rho)
%Simulates pressure and velocity due to a plane wave by means of a
% spherical harmonic synthesis
%
% The sound field is always simulated at 32nd order.
% The incidence direction is given by pw_incidence_azi_deg and 
% pw_incidence_ele_deg.
%
% grid_shape and normal_vector are only required if velocity is to be
% computed.

if nargout == 2
    compute_velocity = 1;
    % assure unit length
    normal_vector = normal_vector ./ vecnorm(normal_vector);
else
    compute_velocity = 0;
    grid_shape = [];
    normal_vector = [];
end

azi_inc_pw = pw_incidence_azi_deg/180*pi;
col_inc_pw = pi/2 - pw_incidence_ele_deg/180*pi;

N = 32; % TODO: take grid size into account

f     = linspace(0, fs/2, taps_pw(1)/2+1).';
omega = 2*pi*f;
k     = omega/c;

% convert to spherical coordinates
[azi_grid, ele_grid, r_grid] = cart2sph(sampling_points(1, :), sampling_points(2, :), sampling_points(3, :));
col_grid = pi/2 - ele_grid;

% ------------ simulate a plane wave using an SH decomposition ------------

% size:                 bins    x  no. of sampl. points x no. of inc. dir.
S              = zeros(length(f), length(azi_grid), length(azi_inc_pw));
S_nm           = zeros(length(f), (N+1)^2, length(azi_inc_pw));

if compute_velocity
   V = zeros(length(f), length(azi_grid), length(azi_inc_pw)); 
end

for index_direction = 1 : length(azi_inc_pw)
    
    for n = 0 : N
       
        bessel = sphbesselj(n, k.*r_grid);
        
        if compute_velocity && strcmp(grid_shape, 'spherical_surface')
            bessel_prime = 1/(2*n+1) * (n * sphbesselj(n-1, k.*r_grid) - (n+1) * sphbesselj(n+1, k.*r_grid));
        end

        for m = -n : n

            % plane wave
            S_nm(:, n^2+n+m+1, index_direction) = sphharm(n, m, col_inc_pw(index_direction), azi_inc_pw(index_direction), 'real');

            % plane wave arrives at coordinate origin at t = 0. Move by taps/2 to middle of buffer. 
            S_nm(:, n^2+n+m+1, index_direction) = S_nm(:, n^2+n+m+1, index_direction) .* exp(-1i .* (0:taps_pw(1)/2).'/taps_pw(1) .* 2*pi .* taps_pw(1)/2);

            S(:, :, index_direction) = S(:, :, index_direction) + S_nm(:, n^2+n+m+1, index_direction) .* 4*pi*1i^n .* bessel .* sphharm(n, m, col_grid, azi_grid, 'real');

            if compute_velocity
                
                if strcmp(grid_shape, 'spherical_surface')
                
                    V(:, :, index_direction) = V(:, :, index_direction) + S_nm(:, n^2+n+m+1, index_direction) .* 4*pi*1i^n .* (-1./(1i*rho*omega)) .* k .* bessel_prime .* sphharm(n, m, col_grid, azi_grid, 'real');
                
                elseif strcmp(grid_shape, 'cubical_surface')
                    
                    V(:, :, index_direction) = V(:, :, index_direction) + S_nm(:, n^2+n+m+1, index_direction) .* 4*pi*1i^n .* normal_vector(1, :) .* (-1./(1i*rho*omega)) .* d_sh_mode(n, m, r_grid, col_grid, azi_grid, k, 'dx');
                    V(:, :, index_direction) = V(:, :, index_direction) + S_nm(:, n^2+n+m+1, index_direction) .* 4*pi*1i^n .* normal_vector(2, :) .* (-1./(1i*rho*omega)) .* d_sh_mode(n, m, r_grid, col_grid, azi_grid, k, 'dy');
                    V(:, :, index_direction) = V(:, :, index_direction) + S_nm(:, n^2+n+m+1, index_direction) .* 4*pi*1i^n .* normal_vector(3, :) .* (-1./(1i*rho*omega)) .* d_sh_mode(n, m, r_grid, col_grid, azi_grid, k, 'dz');
                
                end

            end
            
        end
        
    end
    
end

%fprintf('done.\n\n');

% fix DC
S(1, :, :) = real(S(2, :, :));

% time domain
pressure = ifft(cat(1, S, conj(flipud(S(2:end-1, :, :)))), [], 1, 'symmetric');

% zero pad
pressure = [pressure; zeros(taps_pw(2) - taps_pw(1), size(pressure, 2))];

% do the same for the velocity
if compute_velocity
    
    % fix DC
    V(1, :, :) = real(V(2, :, :));
    velocity = ifft(cat(1, V, conj(flipud(V(2:end-1, :, :)))), [], 1, 'symmetric');
    
    % zero pad
    velocity = [velocity; zeros(taps_pw(2) - taps_pw(1), size(velocity, 2))];

end


end

