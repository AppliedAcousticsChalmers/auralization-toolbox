% (c) by Jens Ahrens, 2024

clear;

addpath('dependencies/');

% sampling grid
grid_file = 'resources/grid_spherical_surface_single_L81.mat';
%grid_file = 'room_data/sound_field_pv_spherical_surface_big_hall_L81.mat'; % comprises both grid and room data (the script disregards the room data)
%grid_file = 'resources/grid_spherical_surface_double_L81.mat';

head_orientation_azimuth_deg = 0; % azimuth in deg, measured counterclockwise

fs = 48000; % sampling frequency

% length of quadrature matrix in time domain
taps_c_nm = 4096; % 4096, good for grids of head size (longer is great, too)

% length of fade-in/fade-out of quadrature matrix
fade_length = 256;

% -------------------------------------------------------------------------

c         = 343; % m/s, speed of sound
rho       = 1.2; % kg/m^3, mass density of air
precision = 'double'; % 'single' (32-bit floating point) or 'double' (64-bit floating point)

% -------------------------------------------------------------------------

% determine parameters for the sample plane wave to be auralized
taps_pw = [1024, 1024+taps_c_nm]; % taps_pw(1): effective length of pw simulation, taps_pw(2): length to zero pad to for the auralization preview

fprintf('\n');

% load the sampling grid (incl. grid_shape and layer_type)
fprintf('Loading sampling grid from file ''%s''.\n\n', grid_file);

load(grid_file);

% dynamic range of the singular values
% (0 < f < 100 Hz, 100 Hz < f < 10 kHz, 10 kHz < f)
if strcmp(grid_shape, 'cubical_volume')
    dynamic_range_dB = [50 40 10];
elseif strcmp(grid_shape, 'cubical_surface')
    dynamic_range_dB = [20 40 10];
elseif strcmp(grid_shape, 'spherical_surface')
    dynamic_range_dB = [20 60 30];
else
    error('Unknown grid_shape.')
end

fprintf('Computing direct auralization matrix for ''%s'' grid and head orientation %d deg. \n', grid_shape, round(head_orientation_azimuth_deg));
fprintf('Regularization: %d | %d | %d dB (f < 100 Hz | f = 100 - 10 kHz) | 10 kHz < f)\n\n', round(dynamic_range_dB));

data_conversion_function = str2func(precision);
    
assert(taps_pw(1) <= taps_c_nm);
assert(rem(taps_c_nm, 2) == 0); % enforce even length

% ----------------------- complement grid data ----------------------------

if strcmp(grid_shape, 'cubical_volume')
    
    % avoid syntax error
    normal_vector = [];

elseif contains(grid_shape, 'surface')

    if strcmp(layer_type, 'double')
  
        % surface normal
        normal_vector = sampling_points_outer - sampling_points_inner;
        normal_vector = normal_vector ./ vecnorm(normal_vector);

        % for file storage and to simplify the syntax
        sampling_points = (sampling_points_inner + sampling_points_outer)/2; 
    
    else % if single-layer surface
        
        if strcmp(grid_shape, 'spherical_surface')   
            % surface normal
            normal_vector = sampling_points/norm(sampling_points(:, 1));
        end

    end

else
    error('Unknown grid shape.');  

end 


%if exist('velocity', 'var')
if strcmp(grid_shape, 'cubical_volume')
    fprintf('Using volumetric pressure to perform the auralization.\n\n');
else
    if strcmp(layer_type, 'single') 
        fprintf('Using surface pressure and particle velocity to perform the auralization.\n\n');
    else
        fprintf('Using double-layer pressure to perform the auralization.\n\n');
    end
end

% ------- create grid of plane wave indicence directions for LS fit -------

% make sure to have more than spatial nodes and that is a square number
no_of_pws = (ceil(sqrt(size(sampling_points, 2))) + 3)^2;

grid_pw = points_on_sphere(no_of_pws);

% convert the pw grid to spherical coordinates
[azi_fliege_rad, ele_fliege_rad, ~] = cart2sph(grid_pw(1, :), grid_pw(2, :), grid_pw(3, :));
azi_fliege_deg = azi_fliege_rad/pi*180;
ele_fliege_deg = ele_fliege_rad/pi*180;

% ------------------------ compute c_l and c_r ----------------------------

% --- eMagLS2 transition frequency ---
no_of_nodes = size(sampling_points, 2);

if strcmp(grid_shape, 'cubical_volume')
    f_transition = 500 * (sqrt(no_of_nodes/2)-1); % approx. 500 Hz * N
elseif strcmp(grid_shape, 'cubical_surface')
    f_transition = 300 * (sqrt(no_of_nodes)-1); % approx. 300 Hz * N
elseif strcmp(grid_shape, 'spherical_surface')
    f_transition = 400 * (sqrt(no_of_nodes)-1); % approx. 500 Hz * N
else 
    error('Unknown grid.');
end

% ------------------ compute auralization matrix --------------------------

fprintf('Loading HRTFs ... ');

% make sure that we have HRTFs for the plane wave incidence directions
hrir_path = 'hrtfs/HRIR_L2702.sofa';
download_hrtfs(hrir_path);

SOFAstart;

hrirs_sofa = SOFAload(hrir_path);

fprintf('done.\n\n');

% TODO: take 'precision' into account
if strcmp(layer_type, 'single') || strcmp(grid_shape, 'cubical_volume')
    [c_l, c_r] = compute_c_direct(head_orientation_azimuth_deg, hrirs_sofa, azi_fliege_deg, ele_fliege_deg, c, taps_c_nm, taps_pw(1), f_transition, fs, dynamic_range_dB, grid_shape, layer_type, normal_vector, rho, sampling_points);
else
    [c_l, c_r] = compute_c_direct(head_orientation_azimuth_deg, hrirs_sofa, azi_fliege_deg, ele_fliege_deg, c, taps_c_nm, taps_pw(1), f_transition, fs, dynamic_range_dB, grid_shape, layer_type, normal_vector, rho, sampling_points, sampling_points_inner, sampling_points_outer);
end

% window the irs
win      = hann(2*fade_length);
fade_in  = win(1:end/2);
fade_out = win(end/2+1:end);

c_l(1:length(fade_in), :, :) = c_l(1:length(fade_in), :, :) .* repmat(fade_in, 1, size(c_l, 2), size(c_l, 3));
c_l(end-length(fade_out)+1:end, :, :) = c_l(end-length(fade_out)+1:end, :, :) .* repmat(fade_out, 1, size(c_l, 2), size(c_l, 3));

c_r(1:length(fade_in), :, :) = c_r(1:length(fade_in), :, :) .* repmat(fade_in, 1, size(c_r, 2), size(c_r, 3));
c_r(end-length(fade_out)+1:end, :, :) = c_r(end-length(fade_out)+1:end, :, :) .* repmat(fade_out, 1, size(c_r, 2), size(c_r, 3));


% ----------------- get sample sound fields for the evaluation ------------
if exist('sampling_points_outer', 'var')
    [~, sampled_sound_field_0 ] = compute_sample_sound_field_for_eq(0,  0, 0, fs, taps_pw, [], grid_shape, normal_vector, c, rho, sampling_points_inner, sampling_points_outer); % sound incidence from straight ahead
    [~, sampled_sound_field_90] = compute_sample_sound_field_for_eq(0, 90, 0, fs, taps_pw, [], grid_shape, normal_vector, c, rho, sampling_points_inner, sampling_points_outer); % sound incidence from the left
else
    [~, sampled_sound_field_0 ] = compute_sample_sound_field_for_eq(0,  0, 0, fs, taps_pw, [], grid_shape, normal_vector, c, rho, sampling_points); % sound incidence from straight ahead
    [~, sampled_sound_field_90] = compute_sample_sound_field_for_eq(0, 90, 0, fs, taps_pw, [], grid_shape, normal_vector, c, rho, sampling_points); % sound incidence from the left
end

% ----- verify auralization of anechoic data against the ground truth -----

fprintf('Computing anechoic auralization data for verification ... ');

sampled_sound_field_0  = double(sampled_sound_field_0 ); % fftfilt requires this
sampled_sound_field_90 = double(sampled_sound_field_90); % fftfilt requires this

brirs_0  = [sum(fftfilt(c_l, sampled_sound_field_0),  2), sum(fftfilt(c_r, sampled_sound_field_0),  2)];
brirs_90 = [sum(fftfilt(c_l, sampled_sound_field_90), 2), sum(fftfilt(c_r, sampled_sound_field_90), 2)];

plot_brirs;

fprintf('done.\n\n');

% auralization preview
create_anechoic_binaural_signals;

% ---------------------------- store everything ---------------------------

if strcmp(grid_shape, 'cubical_volume')
    data_type_string = 'p';
elseif strcmp(layer_type, 'double')
    data_type_string = 'pp';
elseif strcmp(layer_type, 'single')
    data_type_string = 'pv';
else
    error('Something is wrong here.');
end

output_file_name = sprintf('auralization_matrices/auralization_matrix_direct_%s_%s_L%d.mat', data_type_string, grid_shape, size(sampling_points, 2));

fprintf('Storing the auralization matrix in file ''%s'' ... ', output_file_name);

if exist('sampling_points_outer', 'var')
    save(output_file_name, 'c_l', 'c_r', 'fs', 'sampling_points_inner', 'sampling_points_outer', '-v7.3');
else
    save(output_file_name, 'c_l', 'c_r', 'fs', 'sampling_points', '-v7.3');
end

fprintf('done.\n\n');
