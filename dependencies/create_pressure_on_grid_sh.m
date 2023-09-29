function [irs, S_nm] = create_pressure_on_grid_sh(x, y, z, pw_incidence_azi_deg, pw_incidence_ele_deg, c, taps, fs)
%Simulates the pressure due to a plane wave by means of a spherical harmonic
% synthesis
%
% The sound field is always simulated at 32nd order.
% The incidence direction is given by pw_incidence_azi_dega and 
% pw_incidence_ele_deg.

azi_inc_pw  = pw_incidence_azi_deg/180*pi;
col_inc_pw = pi/2 - pw_incidence_ele_deg/180*pi;

N = 32; % TODO: take grid size into account

f = linspace(0, fs/2, taps/2+1).';
k = 2*pi*f/c;

% convert to spherical coordinates
[azi_grid, ele_grid, r_grid] = cart2sph(x, y, z);
col_grid = pi/2 - ele_grid;

% ------------ simulate a plane wave using an SH decomposition ------------

% size:                 bins    x  no. of sampl. points x no. of inc. dir.
S              = zeros(length(f), length(azi_grid), length(azi_inc_pw));
S_nm           = zeros(length(f), (N+1)^2, length(azi_inc_pw));

%fprintf('Computing the sound pressure for direction (azi, ele) = (%d, %d) ... ', round(pw_incidence_azi_deg), round(pw_incidence_ele_deg));

for index_direction = 1 : length(azi_inc_pw)
    
%     if numel(pw_incidence_azi_deg) > 1
%         display_progress(index_direction/length(azi_inc_pw));
%     end
    
    for n = 0 : N
        
%         if numel(pw_incidence_azi_deg) == 1
%             display_progress(n/N);
%         end

        bessel_term = sphbesselj(n, k.*r_grid);

        for m = -n : n

            % plane wave
            S_nm(:, n^2+n+m+1, index_direction) = sphharm(n, m, col_inc_pw(index_direction), azi_inc_pw(index_direction), 'real');

            % plane wave arrives at coordinate origin at t = 0. Move by taps/2 to middle of buffer. 
            S_nm(:, n^2+n+m+1, index_direction) = S_nm(:, n^2+n+m+1) .* exp(-1i .* (0:taps/2).'/taps .* 2*pi .* taps/2);

            S(:, :, index_direction) = S(:, :, index_direction) + S_nm(:, n^2+n+m+1, index_direction) .* 4*pi .* 1i^n .* bessel_term .* sphharm(n, m, col_grid, azi_grid, 'real');

        end
    end
    
end

%fprintf('done.\n\n');

% fix numerical issues
S(1, :, :) = real(S(2, :, :));

% time domain
irs = ifft(cat(1, S, conj(flipud(S(2:end-1, :, :)))), [], 1, 'symmetric');

% % fix data dimensions etc.
% if numel(pw_incidence_azi_deg) == 1
%     grid_data.irs              = squeeze(grid_data.irs);
%     grid_data.S_nm             = squeeze(grid_data.S_nm);
% else
%     grid_data.azi_inc_calibration = azi_inc_pw;
%     grid_data.col_inc_calibration = col_inc_pw;
% end
%
%     if store_output
%         fprintf('Storing output ... ');
%         %save('resources/default_data_volumetric.mat', 'fs', 'irs_verification', 'S_nm', 'x', 'y', 'z', 'c', '-v7.3');
%         save(sprintf('resources/sampling_grids/grid_volumetric_%.2fm_%d_%03ddeg.mat', D, L, azi_deg), 'fs', 'grid_data', 'c', '-v7.3');
%         %save(sprintf('resources/sampling_grids/data_volumetric_%03ddeg_highdens.mat', azi_deg), 'fs', 'irs_verification', 'S_nm', 'x', 'y', 'z', 'c', '-v7.3');
%         fprintf('done.\n\n');
%     else
%         fprintf('\nOutput NOT stored.\n\n');
%     end

