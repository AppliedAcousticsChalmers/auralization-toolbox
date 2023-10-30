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

N    = 8; % (N+1)^2 directions of incidence are used in the computation
taps = 1024; % good for grids of head size

% regularization
reg_parameter = [10^(20/60) 10^(20/60)]; % gentle: 20 dB SNR (sigma_n^2 / sigma_s^2)
fc_reg = [50 16000]; % frequency ranges in which to regularize

c         = 343; % m/s, speed of sound
rho       = 1.2; % kg/m^3, mass density of air
precision = 'single'; % 'single' (32-bit flating point) or 'double' (64-bit floating point)

% -------------------------------------------------------------------------

fprintf('\n');
fprintf('Processing ''%s'' data of room ''%s''.\n\n', grid_shape, room);

data_conversion_function = str2func(precision);
    
assert(rem(taps, 2) == 0); % enforce even length

% ------------------------------ get HRTFs --------------------------------

% create plane wave data for least-squares fit

% get a Fliege grid of sample incidence directions
grid_pw = load(sprintf('resources/fliege_grids/%d.txt', (N+1)^2));

[azi_pw_deg, ele_pw_deg, ~] = cart2sph(grid_pw(:, 1).', grid_pw(:, 2).', grid_pw(:, 3).');
azi_pw_deg = azi_pw_deg/pi*180;
ele_pw_deg = ele_pw_deg/pi*180;

% make sure that we have HRTFs for the desired angles
hrir_path = 'hrtfs/HRIR_L2702.sofa';
download_hrtfs(hrir_path);

SOFAstart;
hrirs_sofa = SOFAload(hrir_path);

% rotate HRTF set to account for the head orientation
hrirs_sofa.SourcePosition(:, 1) = hrirs_sofa.SourcePosition(:, 1) + head_orientation_deg;

[azi_pw_deg, ele_pw_deg, indices_hrirs_calibration] = pick_angles_according_to_hrtfs(hrirs_sofa, azi_pw_deg, ele_pw_deg);

% ---------------------- compute calibration data -------------------------

if strcmp(grid_shape, 'cubical_volume')
    
    % make sure that we're getting an overdetermined equation system
    fprintf('The minimum required number of sampling points for the desired N is %d. %d are available.\n\n', (N+1)^2, size(sampling_points, 2));

    % dimensions:           taps  x  no. of sampling points  x  no. of plane waves
    sampled_sound_field_tmp = zeros(taps, size(sampling_points, 2), length(azi_pw_deg));
 
    display_progress('Computing the calibration data:');
    
    for index = 1 : length(azi_pw_deg)

        display_progress(index/length(azi_pw_deg));

        sampled_sound_field_tmp(:, :, index) = compute_plane_wave_on_grid(sampling_points(1, :), sampling_points(2, :), sampling_points(3, :), azi_pw_deg(index), ele_pw_deg(index), c, taps, fs);

    end
    
elseif strcmp(grid_shape, 'spherical_surface') || strcmp(grid_shape, 'cubical_surface')
    
    if exist('velocity', 'var')
                
        fprintf('Using pressure and particle velocity to perform the auralization.\n\n');
        
        % make sure that we're getting an overdetermined equation system
        fprintf('The minimum required number of pairs of sampling points for the desired N is %d. %d are available.\n\n', (N+1)^2, size(sampling_points, 2));
    
        % dimensions:           taps  x  no. of sampling points  x  no. of plane waves
        sampled_sound_field_tmp = zeros(taps, size(sampling_points, 2), length(azi_pw_deg));
    
    else
        
        fprintf('Using double pressure layer to perform the auralization.\n\n');
        
        % make sure that we're getting an overdetermined equation system
        fprintf('The minimum required number of pairs of sampling points for the desired N is %d. %d are available.\n\n', (N+1)^2, size(sampling_points_outer, 2));
    
        % dimensions:           taps  x  no. of sampling points  x  no. of plane waves
        sampled_sound_field_tmp = zeros(taps, size(sampling_points_outer, 2), length(azi_pw_deg));
    
    end
    
    display_progress('Computing the calibration data:');
    
    for index = 1 : length(azi_pw_deg)
        
        display_progress(index/length(azi_pw_deg));
        
        if exist('velocity', 'var')
            
            [pressure_tmp, velocity_tmp] = compute_plane_wave_on_grid(sampling_points(1, :), sampling_points(2, :), sampling_points(3, :), azi_pw_deg(index), ele_pw_deg(index), c, taps, fs, grid_shape, normal_vector, rho);
            
            % compute the signals captured by the virtual cardioid transducers
            sampled_sound_field_tmp(:, :, index) = pressure_tmp - rho*c .* velocity_tmp;
            
        else
            
            pressure_outer_tmp = compute_plane_wave_on_grid(sampling_points_outer(1, :), sampling_points_outer(2, :), sampling_points_outer(3, :), azi_pw_deg(index), ele_pw_deg(index), c, taps, fs);
            pressure_inner_tmp = compute_plane_wave_on_grid(sampling_points_inner(1, :), sampling_points_inner(2, :), sampling_points_inner(3, :), azi_pw_deg(index), ele_pw_deg(index), c, taps, fs);
       
            % TODO: Implement computation of particle velocity
            [~, pressure_tmp] = compute_cardioid_from_pressure(sampling_points_inner, sampling_points_outer, pressure_inner_tmp, pressure_outer_tmp, fs, c);
    
            sampled_sound_field_tmp(:, :, index) = pressure_tmp;
            
        end
    end
    
else
    error('Unknown grid shape.');  
end

fprintf('\n\n');

% ------------------- compute quadrature matrix C_nm ----------------------

fprintf('Computing the quadrature matrices ... ');

% the dimensions of C_nm are (no. of frequency bins x no. of sampling positions x XXXXXXXXXXXX)        
[C_nm_left, condition_number] = get_c_nm_direct(sampled_sound_field_tmp, squeeze(hrirs_sofa.Data.IR(indices_hrirs_calibration, 1, :)).', fs, fc_reg, reg_parameter);
C_nm_right = get_c_nm_direct(sampled_sound_field_tmp, squeeze(hrirs_sofa.Data.IR(indices_hrirs_calibration, 2, :)).', fs, fc_reg, reg_parameter);

fprintf('done. The maximum condition number is %d.\n\n', round(max(condition_number)));

% conversion to time domain
c_nm_left  = ifft(cat(1, C_nm_left,  flipud(conj(C_nm_left( 2:end-1, :, :)))), [], 'symmetric');
c_nm_right = ifft(cat(1, C_nm_right, flipud(conj(C_nm_right(2:end-1, :, :)))), [], 'symmetric');

% enforce causality (causes a time delay of taps/2)
c_nm_left  = circshift(c_nm_left,  size(c_nm_left,  1)/2, 1);
c_nm_right = circshift(c_nm_right, size(c_nm_right, 1)/2, 1);

% -------------------------- binaural rendering ---------------------------

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

sampled_sound_field = double(sampled_sound_field); % fftfilt requires this

brirs = [sum(fftfilt(c_nm_left, sampled_sound_field), 2), sum(fftfilt(c_nm_right, sampled_sound_field), 2)];

% get an example audio signal
[sig, fs_sig] = audioread('resources/drum_loop_48k.wav');

assert(fs == fs_sig);

% render audio signal binaurally
out_binaural = fftfilt(brirs, sig);

% normalize
out_binaural = out_binaural ./ max(abs(out_binaural(:)));

% store binaural signal
audiowrite('out_binaural_direct.wav', out_binaural, fs);



