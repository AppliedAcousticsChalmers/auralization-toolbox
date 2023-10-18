% This script convert surface sampled data to an ambisonic
% representation which in turn is rendered binaurally.
%
% (c) by Jens Ahrens, 2023

clear;

addpath('dependencies/');

% load the variables fs, sampling_points, grid_shape, pressure, velocity, room
load('resources/sound_field_p_cubical_volume_living_room.mat');
%load('resources/sound_field_p_cubical_volume_big_hall.mat');
%load('resources/sound_field_pp_spherical_surface_living_room.mat');
%load('resources/sound_field_pp_spherical_surface_big_hall.mat');
%load('resources/sound_field_pp_cubical_surface_living_room.mat');
%load('resources/sound_field_pp_cubical_surface_big_hall.mat');
%load('resources/sound_field_pv_spherical_surface_living_room.mat');
%load('resources/sound_field_pv_spherical_surface_big_hall.mat');
%load('resources/sound_field_pv_cubical_surface_living_room.mat');
%load('resources/sound_field_pv_cubical_surface_big_hall.mat');

head_orientation_deg = 0; % azimuth in deg, counterclockwise

N    = 8; % order of the intermediate spherical harmonic representation
taps = 1024; % good for grids of head size

% Regularization of the SH matrix. No element may be smaller in magnitude
% than -50 dB below the maximum.
dynamic_range_dB = 50; % dB

c         = 343; % m/s, speed of sound
rho       = 1.2; % kg/m^3, mass density of air
precision = 'single'; % 'single' (32-bit flating point) or 'double' (64-bit floating point)

% -------------------------------------------------------------------------

fprintf('\n');
fprintf('Processing ''%s'' data of room ''%s''.\n\n', grid_shape, room);

data_conversion_function = str2func(precision);
    
assert(rem(taps, 2) == 0); % enforce even length

% ------------------- compute quadrature matrix C_nm ---------------------

if strcmp(grid_shape, 'cubical_volume')
    
    % make sure that we're getting an overdetermined equation system
    fprintf('The minimum required number of sampling points for the desired N is %d. %d are available.\n\n', (N+1)^2, size(sampling_points, 2));

    % the dimensions of C_nm are (no. of frequency bins x no. of sampling positions x (N+1)^2)
    [C_nm, condition_number] = get_c_nm_volumetric(taps, sampling_points, fs, N, c, dynamic_range_dB, 'real', precision);

else % surface sampling
    
    % single layer
    if exist('velocity', 'var')
        
        fprintf('Using pressure and particle velocity to perform the SH decompostion.\n\n');
    
        % make sure that we're getting an overdetermined equation system
        fprintf('The minimum required number of pairs of sampling points for the desired N is %d. %d are available.\n\n', (N+1)^2, size(sampling_points, 2));
    
        [azi_pressure, ele_pressure, r_pressure] = cart2sph(sampling_points(1, :).', sampling_points(2, :).', sampling_points(3, :).');
        col_pressure = pi/2 - ele_pressure;
        
        r_velocity = r_pressure;
        
    % double layer
    else 
        
        fprintf('Using double pressure layer to perform the SH decompostion.\n\n');
        
        % make sure that we're getting an overdetermined equation system
        fprintf('The minimum required number of pairs of sampling points for the desired N is %d. %d are available.\n\n', (N+1)^2, size(sampling_points_outer, 2));
    
        % for measuring the pressure
        [azi_pressure, ele_pressure, r_pressure] = cart2sph(sampling_points_outer(1, :).', sampling_points_outer(2, :).', sampling_points_outer(3, :).');
        col_pressure = pi/2 - ele_pressure;

        % for measuring the velocity
        [~, ~, r_velocity] = cart2sph(sampling_points_inner(1, :).', sampling_points_inner(2, :).', sampling_points_inner(3, :).'); 
        r_velocity = (r_pressure + r_velocity)/2; % mid-point between layers
        
        if strcmp(grid_shape, 'cubical_surface') 
        
            % compute xyz compontents of the normal vector in direction of the derivative
            normal_vector = sampling_points_outer - sampling_points_inner;
            normal_vector = normal_vector ./ vecnorm(normal_vector);
        
        end
        
    end
    
    % finally, compute C_nm
    if strcmp(grid_shape, 'spherical_surface')
        [C_nm, condition_number] = get_c_nm_surface_radial(taps, r_pressure(1), r_velocity(1), azi_pressure, col_pressure, fs, N, c, dynamic_range_dB, precision); 
    elseif strcmp(grid_shape, 'cubical_surface') 
        [C_nm, condition_number] = get_c_nm_surface_cartesian(taps, r_pressure, r_velocity, normal_vector, azi_pressure, col_pressure, fs, N, c, dynamic_range_dB, precision);
    else
        error('Unknown grid shape.');
    end
    
end % if cubical or spherical sampling
    
fprintf('The maximum condition number is %d.\n\n', round(max(condition_number)));

% conversion to time domain
c_nm = ifft(cat(1, C_nm, flipud(conj(C_nm(2:end-1, :, :)))), [], 'symmetric');

% enforce causality (causes a time delay of taps/2)
c_nm = circshift(c_nm, size(c_nm, 1)/2, 1);

% ------------------------- SH decomposition ------------------------------

% prepare the sound field data for SH decomposition
if strcmp(grid_shape, 'spherical_surface') || strcmp(grid_shape, 'cubical_surface')
    
    if exist('velocity', 'var')
        % compute the signals captured by the virtual cardioid transducers
        sampled_sound_field = pressure - rho*c .* velocity;
    else
        % compute the signals captured by the virtual cardioid transducers
        [~, sampled_sound_field] = compute_cardioid_from_pressure(sampling_points_inner, sampling_points_outer, pressure_inner, pressure_outer, fs, c);
    end
    
else
    
    sampled_sound_field = pressure;
    
end

% ambisonic representation of the room data 
s_nm = get_sh_coefficients_t(sampled_sound_field, c_nm, N);

% % get an example audio signal
% [sig, fs_sig] = audioread('resources/drum_loop_48k.wav');
% 
% assert(fs == fs_sig);
% 
% % create a running ambisonic signal from the ambisonic room impulse response
% out_ambisonic = fftfilt(s_nm, sig);
% 
% % normalize
% out_ambisonic = out_ambisonic ./ max(abs(out_ambisonic(:)));
% 
% % store binaural signal
% audiowrite('out_ambisonic.wav', out_ambisonic, fs);

% -------------------------- binaural rendering ---------------------------

% HRIRs are loaded inside the function
brirs = render_ambi_signals_binaurally_t(s_nm, head_orientation_deg/180*pi, N, 'transform_integral');

% get an example audio signal
[sig, fs_sig] = audioread('resources/drum_loop_48k.wav');

assert(fs == fs_sig);

% render audio signal binaurally
out_binaural = fftfilt(brirs, sig);

% normalize
out_binaural = out_binaural ./ max(abs(out_binaural(:)));

% store binaural signal
audiowrite('out_binaural_via_ambisonics.wav', out_binaural, fs);


