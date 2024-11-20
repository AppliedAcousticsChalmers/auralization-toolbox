clear;

addpath('dependencies/');

% This script works for both ambisonic and direct auralization.

% auralization_matrix_file (comprises also the sampling grid coordinates)
auralization_matrix_file = 'auralization_matrices/auralization_matrix_ambisonics_pv_spherical_surface_L81.mat';

% sound field data
simulation_data_file     = 'room_data/sound_field_pv_spherical_surface_big_hall_L81.mat'; 

audio_file               = 'resources/drum_loop_48k.wav';

% head_orientation_azimuth_deg is only taken into account with ambisonic
% auralization; the head orientation is baked into the auralization matrix
% in direct auralization
head_orientation_azimuth_deg = 0;

% -------------------------------------------------------------------------

c   = 343; % m/s, speed of sound
rho = 1.2; % kg/m^3, mass density of air

auralization_matrix = load(auralization_matrix_file);

fprintf('\n');

if exist('simulation_data_file', 'var')
    
    simulation_data = load(simulation_data_file);
    
    assert(simulation_data.fs == auralization_matrix.fs);
    
    fprintf('Auralizing the data in file ''%s''\n', simulation_data_file);

else 

    fprintf('Auralizing anechoic data\n');

    % avoid syntax errors
    simulation_data.pressure = 0;

end

fprintf('with the auralization matrix in file ''%s''.\n\n', auralization_matrix_file);

% -------------------- figure out the grid geometry -----------------------

fprintf('Preparing the simulation data ... ');

% single layer surface
if isfield(simulation_data, 'velocity')
    sampled_sound_field = simulation_data.pressure - rho*c .* simulation_data.velocity;
% double layer
elseif isfield(simulation_data, 'pressure_inner')
    [~, sampled_sound_field] = compute_cardioid_from_pressure(simulation_data.sampling_points_inner, simulation_data.sampling_points_outer, simulation_data.pressure_inner, simulation_data.pressure_outer, fs, c);
% volumetric grid
else
    sampled_sound_field = simulation_data.pressure;
end

fprintf('done.\n\n');

% --- if ambisonic auralization ---
if isfield(auralization_matrix, 'c_nm')

    % determine ambisonic order
    N = sqrt(size(auralization_matrix.c_nm, 3))-1;

    fprintf('Detected ambisonic auralization matrix of order %d.\n\n', N);
    fprintf('Computing the binaural signal ... ');

    s_nm  = get_sh_coefficients_t(sampled_sound_field, auralization_matrix.c_nm, N);

    % load MagLS-equalized HRTFs
    hrirs = load(sprintf('hrtfs/hrirs_ku100_magls_sh_N%d.mat', N));

    % finally, auralize
    brirs = render_ambi_signals_binaurally_t(s_nm, hrirs, auralization_matrix.eq_ir, head_orientation_azimuth_deg/180*pi, N);

% --- direct auralization ---
elseif isfield(auralization_matrix, 'c_l')

    fprintf('Detected direct auralization matrix. head_orientation_azimuth_deg is not taken into account.\n\n');
    fprintf('Computing the binaural signal ... ');

    % fftfilt requires this
    sampled_sound_field     = double(sampled_sound_field); 
    auralization_matrix.c_l = double(auralization_matrix.c_l);
    auralization_matrix.c_r = double(auralization_matrix.c_r);

    % apply the auralization matrix
    brirs = [sum(fftfilt(auralization_matrix.c_l, sampled_sound_field),  2), sum(fftfilt(auralization_matrix.c_r, sampled_sound_field),  2)];

else

    error('Cannot identify the auralization matrix.');

end

fprintf('done.\n\n'); % 'Computing the binaural signal ...'

% ----------------------- create audio example ----------------------------

[sig, fs_sig] = audioread(audio_file);

assert(simulation_data.fs == fs_sig);

% render audio signal binaurally
auralized_binaural = fftfilt(brirs, sig(:, 1));

% normalize
auralized_binaural = auralized_binaural ./ max(abs(auralized_binaural(:))) * .99;

fprintf('Storing binaural audio example in file ''auralization_binaural.wav'' ... ');

audiowrite('auralization_binaural.wav', auralized_binaural, fs_sig);

fprintf('done.\n\n');

% --- if ambisonic auralization ---
if isfield(auralization_matrix, 'c_nm')
    
    fprintf('Creating ambisonic audio example ... ');

    % create ambisonic signal
    auralized_ambisonic = fftfilt(s_nm, sig(:, 1));

    % normalize
    auralized_ambisonic = auralized_ambisonic ./ max(abs(auralized_ambisonic(:))) * .99;

    fprintf('done.\n\n');

    fprintf('Storing ambisonic audio example in file ''auralization_ambisonic_ACN_N3D.wav'' ... ');

    audiowrite('auralization_ambisonic_ACN_N3D.wav', auralized_ambisonic, fs_sig);

    fprintf('done.\n\n');
    
end



