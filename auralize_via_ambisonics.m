% This script convert surface sampled data to an ambisonic
% representation which in turn is rendered binaurally.
%
% (c) by Jens Ahrens, 2023

clear;

addpath('dependencies/');

% load the variables fs, grid_data, grid_shape, irs, room
%load('resources/sound_field_cubical_volume_living_room.mat');
load('resources/sound_field_cubical_volume_big_hall.mat');
%load('resources/sound_field_spherical_surface_living_room.mat');
%load('resources/sound_field_spherical_surface_big_hall.mat');
%load('resources/sound_field_cubical_surface_living_room.mat');
%load('resources/sound_field_cubical_surface_big_hall.mat');

head_orientation_deg = 45; % azimuth in deg, counterclockwise

N    = 8; % order of the intermediate spherical harmonic representation
taps = 1024; % good for grids of head size

% Regularization of the SH matrix. No element may be smaller in magnitude
% than -50 dB below the maximum.
dynamic_range_dB = 50; % dB

c         = 343; % m/s
precision = 'single'; % 'single' (32-bit flating point) or 'double' (64-bit floating point)

% -------------------------------------------------------------------------

fprintf('\n');
fprintf('Processing ''%s'' data of room ''%s''.\n\n', grid_shape, room);

data_conversion_function = str2func(precision);
    
assert(rem(taps, 2) == 0); % enforce even length

% ------------------- compute quadrature matrix C_nm ---------------------

if strcmp(grid_shape, 'cubical_volume')
    
    % make sure that we're getting an overdetermined equation system
    fprintf('The minimum required number of sampling points for the desired N is %d. You chose %d.\n\n', (N+1)^2, size(grid_data.sampling_points, 2));

    % the dimensions of C_nm are (no. of frequency bins x no. of sampling positions x (N+1)^2)
    [C_nm, condition_number] = get_c_nm_volumetric(taps, grid_data, fs, N, c, dynamic_range_dB, 'real', precision);

elseif strcmp(grid_shape, 'spherical_surface')
    
    [C_nm, condition_number] = get_c_nm_surface_radial(taps, grid_data, fs, N, c, dynamic_range_dB, precision); 

elseif strcmp(grid_shape, 'cubical_surface')

    [C_nm, condition_number] = get_c_nm_surface_cartesian(taps, grid_data, fs, N, c, dynamic_range_dB, precision);

else
    error('Unknown grid shape.');
end
    
fprintf('The maximum condition number is %d.\n\n', round(max(condition_number)));

% conversion to time domain
c_nm = ifft(cat(1, C_nm, flipud(conj(C_nm(2:end-1, :, :)))), [], 'symmetric');

% enforce causality (causes a time delay of taps/2)
c_nm = circshift(c_nm, size(c_nm, 1)/2, 1);

% ------------------------- SH decomposition ------------------------------

if strcmp(grid_shape, 'spherical_surface') || strcmp(grid_shape, 'cubical_surface')
    
    % compute the signals captured by the virtual cardioid transducers
    [~, irs] = compute_cardioid_from_pressure(grid_data, irs_outer, irs_inner, fs, c);
    
end

% ambisonic representation of the room data
s_nm = get_sh_coefficients_t(irs, c_nm, N);

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


