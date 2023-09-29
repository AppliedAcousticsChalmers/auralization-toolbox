% This script convert surface sampled data to an ambisonic
% representation which in turn is rendered binaurally.
%
% (c) by Jens Ahrens, 2023

clear;

addpath('dependencies/');

% load the variables fs, grid_data, grid_shape, irs, room
load('resources/sound_field_cubical_volume_living_room.mat');
%load('resources/sound_field_cubical_volume_big_hall.mat');
%load('resources/sound_field_spherical_surface_living_room.mat');
%load('resources/sound_field_spherical_surface_big_hall.mat');
%load('resources/sound_field_cubical_surface_living_room.mat');
%load('resources/sound_field_cubical_surface_big_hall.mat');

head_orientation_deg = 45; % azimuth in deg

N    = 8; % (N+1)^2 directions of incidence are used in the computation
taps = 1024; % good for grids of head size

% regularization
reg_parameter = [10^(-20/60) 10^(-20/60)]; % gentle: 20 dB SNR (sigma_n^2 / sigma_s^2)
fc_reg = [50 16000]; % frequency ranges is which to regularize

c         = 343; % m/s
precision = 'single'; % 'single' (32-bit flating point) or 'double' (64-bit floating point)

% -------------------------------------------------------------------------

fprintf('\n');
fprintf('Processing ''%s'' data of room ''%s''.\n\n', grid_shape, room);

data_conversion_function = str2func(precision);
    
assert(rem(taps, 2) == 0); % enforce even length

% ------------------- compute quadrature matrix C_nm ----------------------

% create plane wave data for least-squares fit

% get a Fliege grid of sample incidence directions
grid_pw = load(sprintf('resources/fliege_grids/%d.txt', (N+1)^2));

[azi_pw_deg, ele_pw_deg, ~] = cart2sph(grid_pw(:, 1).', grid_pw(:, 2).', grid_pw(:, 3).');
azi_pw_deg = azi_pw_deg/pi*180;
ele_pw_deg = ele_pw_deg/pi*180;

% download HRIRs (skipped automatically if the HRIR dataset already exists)
hrir_path = 'hrtfs/HRIR_L2702.sofa';

if ~isfile(hrir_path)
    mkdir hrtfs;
    fprintf('Can''t find the HRTFs. '); 
    fprintf('Downloading them from https://zenodo.org/record/3928297 ... ');
    websave(hrir_path, 'https://zenodo.org/record/3928297/files/HRIR_L2702.sofa?download=1');
    fprintf('done.\n\n');
end

% make sure that we have HRTFs for the desied angles
SOFAstart;
hrirs_sofa = SOFAload('hrtfs/HRIR_L2702.sofa');

% rotate HRTF set to account for the head orientation
hrirs_sofa.SourcePosition(:, 1) = hrirs_sofa.SourcePosition(:, 1) + head_orientation_deg;

[azi_pw_deg, ele_pw_deg, indices_hrirs_calibration] = pick_angles_according_to_hrtfs(hrirs_sofa, azi_pw_deg, ele_pw_deg);

if strcmp(grid_shape, 'cubical_volume')
    
    % make sure that we're getting an overdetermined equation system
    fprintf('The minimum required number of sampling points for the desired N is %d. You chose %d.\n\n', (N+1)^2, size(grid_data.sampling_points, 2));

    % dimensions:           taps  x  no. of sampling points  x  no. of plane waves
    irs_calibration = zeros(taps, size(grid_data.sampling_points, 2), length(azi_pw_deg));
 
    display_progress('Computing the calibration data:');
    
    for index = 1 : length(azi_pw_deg)

        display_progress(index/length(azi_pw_deg));

        irs_calibration(:, :, index) = create_pressure_on_grid_sh(grid_data.sampling_points(1, :), grid_data.sampling_points(2, :), grid_data.sampling_points(3, :), azi_pw_deg(index), ele_pw_deg(index), c, taps, fs);

    end
    
elseif strcmp(grid_shape, 'spherical_surface') || strcmp(grid_shape, 'cubical_surface')
    
    % make sure that we're getting an overdetermined equation system
    fprintf('The minimum required number of sampling points for the desired N is %d. You chose %d.\n\n', (N+1)^2, size(grid_data.sampling_points_outer, 2));

    % dimensions:           taps  x  no. of sampling points  x  no. of plane waves
    irs_calibration = zeros(taps, size(grid_data.sampling_points_outer, 2), length(azi_pw_deg));
    
    display_progress('Computing the calibration data:');
    
    for index = 1 : length(azi_pw_deg)
        
        display_progress(index/length(azi_pw_deg));
        
        irs_outer = create_pressure_on_grid_sh(grid_data.sampling_points_outer(1, :), grid_data.sampling_points_outer(2, :), grid_data.sampling_points_outer(3, :), azi_pw_deg(index), ele_pw_deg(index), c, taps, fs);
        irs_inner = create_pressure_on_grid_sh(grid_data.sampling_points_inner(1, :), grid_data.sampling_points_inner(2, :), grid_data.sampling_points_inner(3, :), azi_pw_deg(index), ele_pw_deg(index), c, taps, fs);
       
        [~, irs_tmp] = compute_cardioid_from_pressure(grid_data, irs_outer, irs_inner, fs, c);
        
        irs_calibration(:, :, index) = irs_tmp;
    end
    
else
    error('Unknown grid shape.');  
end

fprintf('\n\n');

fprintf('Computing the quadrature matrices ... ');

% the dimensions of C_nm are (no. of frequency bins x no. of sampling positions x XXXXXXXXXXXX)        
[C_nm_left, condition_number] = get_c_nm_direct(irs_calibration, squeeze(hrirs_sofa.Data.IR(indices_hrirs_calibration, 1, :)).', fs, fc_reg, reg_parameter);
C_nm_right = get_c_nm_direct(irs_calibration, squeeze(hrirs_sofa.Data.IR(indices_hrirs_calibration, 2, :)).', fs, fc_reg, reg_parameter);

fprintf('done. The maximum condition number is %d.\n\n', round(max(condition_number)));

% conversion to time domain
c_nm_left  = ifft(cat(1, C_nm_left,  flipud(conj(C_nm_left( 2:end-1, :, :)))), [], 'symmetric');
c_nm_right = ifft(cat(1, C_nm_right, flipud(conj(C_nm_right(2:end-1, :, :)))), [], 'symmetric');

% enforce causality (causes a time delay of taps/2)
c_nm_left  = circshift(c_nm_left,  size(c_nm_left,  1)/2, 1);
c_nm_right = circshift(c_nm_right, size(c_nm_right, 1)/2, 1);

% -------------------------- binaural rendering ---------------------------

if strcmp(grid_shape, 'spherical_surface') || strcmp(grid_shape, 'cubical_surface')
    [~, irs] = compute_cardioid_from_pressure(grid_data, irs_outer, irs_inner, fs, c);
end

irs = double(irs); % fftfilt requires this

brirs = [sum(fftfilt(c_nm_left, irs), 2), sum(fftfilt(c_nm_right, irs), 2)];

% get an example audio signal
[sig, fs_sig] = audioread('resources/drum_loop_48k.wav');

assert(fs == fs_sig);

% render audio signal binaurally
out_binaural = fftfilt(brirs, sig);

% normalize
out_binaural = out_binaural ./ max(abs(out_binaural(:)));

% store binaural signal
audiowrite('out_binaural_direct.wav', out_binaural, fs);




