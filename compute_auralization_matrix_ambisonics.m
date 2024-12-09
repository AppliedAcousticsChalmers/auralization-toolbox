% (c) by Jens Ahrens, 2024

clear;

addpath('dependencies/');

% sampling grid
grid_file = 'resources/grid_spherical_surface_L81.mat';
%grid_file = 'room_data/sound_field_pv_spherical_surface_big_hall_L81.mat'; % comprises both grid and room data (the script disregards the room data)

N = 7; % desired spherical harmonic order; make sure to check Tab. 1 in 
       % the JAES manuscript to understand what order a given grid supports

head_orientation_azimuth_deg = 0; % this is only for computing the verificatiom data

fs = 48000; % sampling frequency

% length of quadrature matrix in time domain
taps_c_nm = 1024; % 1024, good for grids of head size (longer is great, too)
taps_pw   = [1024 2048]; % taps_pw(1): effective length of pw simulation, taps_pw(2): length to zero pad to in the evaluation

% use the following with cubical volume grids to avoid the faint echo on
% the contralateral side
%taps_c_nm = 4096;
%taps_pw   = [1024 2*4096];

% -------------------------------------------------------------------------

c         = 343; % m/s, speed of sound
rho       = 1.2; % kg/m^3, mass density of air
precision = 'single'; % 'single' (32-bit floating point) or 'double' (64-bit floating point)

% -------------------------------------------------------------------------

hrtf_type = 'magls';

fprintf('\n');

% load the sampling grid (incl. grid_shape and layer_type)
fprintf('Loading sampling grid from file ''%s''.\n\n', grid_file);

load(grid_file);

fprintf('Computing ambisonic auralization matrix for ''%s'' grid. N = %d.\n\n', grid_shape, N);

% Regularization of the SH matrix, dynamic range of singular values
if (strcmp(grid_shape, 'spherical_surface') && N > 4)
    dynamic_range_dB = 70; 
else
    dynamic_range_dB = 40;
end

fprintf('Limiting the dynamic range of the singular values to %d dB.\n\n', round(dynamic_range_dB));


data_conversion_function = str2func(precision);

assert(taps_pw(1) <= taps_c_nm);
assert(rem(taps_c_nm, 2) == 0); % enforce even length

% ------------------- compute quadrature matrix C_nm ---------------------

if strcmp(grid_shape, 'cubical_volume')
      
    fprintf(2, ['Consider uncommenting the lines 24-25 if you are using a cubical ' ...
                'volumetric grid, and you are unsatisfied with the result.\n\n']);

    % avoid syntax error
    normal_vector = [];

    % check if overdetermined equation system
    check_equation_system((N+1)^2, size(sampling_points, 2), 'volume');

    % the dimensions of C_nm are (no. of frequency bins x no. of sampling positions x (N+1)^2)
    C_nm = get_c_nm_volumetric(taps_c_nm, sampling_points, fs, N, c, dynamic_range_dB, 'real', precision);
    
elseif contains(grid_shape, 'surface')

    if strcmp(grid_shape, 'spherical_surface')
              
        if strcmp(layer_type, 'single') 

            % surface normal
            normal_vector = sampling_points/norm(sampling_points(:, 1));
           
        elseif strcmp(layer_type, 'double')
    
            % surface normal
            normal_vector = sampling_points_outer - sampling_points_inner;
            normal_vector = normal_vector ./ vecnorm(normal_vector);

            % for file storage and to simplify the syntax
            sampling_points = (sampling_points_inner + sampling_points_outer)/2; 

        end

        % check if overdetermined equation system
        check_equation_system((N+1)^2, size(sampling_points, 2), 'surface');
        
    elseif strcmp(grid_shape, 'cubical_surface')
        
        % surface normal
        normal_vector = sampling_points_outer - sampling_points_inner;
        normal_vector = normal_vector ./ vecnorm(normal_vector);
       
        % cubical single-layer surface
        if strcmp(layer_type, 'single') 

            sampling_points = sampling_points_outer;

            clear sampling_points_inner sampling_points_outer;     

        elseif strcmp(layer_type, 'double')
        
            % for file storage and to simplify the syntax
            sampling_points = (sampling_points_inner + sampling_points_outer)/2; 

        end
        
        % check if overdetermined equation system
        check_equation_system((N+1)^2, size(sampling_points, 2), 'surface');

    else 
        error('Unknown grid_shape.');
    end

    % single layer
    if strcmp(layer_type, 'single')
        
        fprintf('Using pressure and particle velocity to perform the SH decompostion.\n\n');
        
        [azi_pressure, ele_pressure, r_pressure] = cart2sph(sampling_points(1, :).', sampling_points(2, :).', sampling_points(3, :).');
        col_pressure = pi/2 - ele_pressure;
        
        azi_velocity = azi_pressure;
        col_velocity = col_pressure;
        r_velocity = r_pressure;
        
    % double layer
    else 
        
        fprintf('Using double pressure layer to perform the SH decompostion.\n\n');
            
        % for measuring the pressure
        [azi_pressure, ele_pressure, r_pressure] = cart2sph(sampling_points_outer(1, :).', sampling_points_outer(2, :).', sampling_points_outer(3, :).');
        col_pressure = pi/2 - ele_pressure;

        % mid-point between layers
        sampling_points_velocity = (sampling_points_inner + sampling_points_outer)/2;

        % for measuring the velocity
        [azi_velocity, ele_velocity, r_velocity] = cart2sph(sampling_points_velocity(1, :).', sampling_points_velocity(2, :).', sampling_points_velocity(3, :).'); 
        col_velocity = pi/2 - ele_velocity; 
                
    end
    
    % finally, compute C_nm
    if strcmp(grid_shape, 'spherical_surface')
        
        % the dimensions of C_nm are (no. of frequency bins x no. of sampling positions x (N+1)^2)
        C_nm = get_c_nm_surface_radial(taps_c_nm, r_pressure(1), r_velocity(1), azi_pressure, col_pressure, fs, N, c, dynamic_range_dB, precision); 

    elseif strcmp(grid_shape, 'cubical_surface') 

        % the dimensions of C_nm are (no. of frequency bins x no. of sampling positions x (N+1)^2)
        C_nm = get_c_nm_surface_cartesian(taps_c_nm, r_pressure, azi_pressure, col_pressure, r_velocity, azi_velocity, col_velocity, normal_vector, fs, N, c, dynamic_range_dB, precision);
      
    else
       error('Unknown grid shape.');
    end

else
    error('Unknown grid shape.');  

end % if cubical or spherical sampling
    
% conversion to time domain
c_nm = ifft(cat(1, C_nm, flipud(conj(C_nm(2:end-1, :, :)))), [], 'symmetric');

% enforce causality (causes a time delay of taps/2)
c_nm = circshift(c_nm, size(c_nm, 1)/2, 1);

% window the irs, just in case
win      = hann(256);
fade_in  = win(1:end/2);
fade_out = win(end/2+1:end);

c_nm(1:length(fade_in), :, :) = c_nm(1:length(fade_in), :, :) .* repmat(fade_in, 1, size(c_nm, 2), size(c_nm, 3));
c_nm(end-length(fade_out)+1:end, :, :) = c_nm(end-length(fade_out)+1:end, :, :) .* repmat(fade_out, 1, size(c_nm, 2), size(c_nm, 3));

% ---------------------------- compute EQ ---------------------------------

fprintf('Computing the sound field EQ ... ');

% compute a plane wave that impinges from straight ahead
if exist('sampling_points_outer', 'var')
    s_nm = compute_sample_sound_field_for_eq(c_nm, 0, 0, fs, taps_pw, N, grid_shape, normal_vector, c, rho, sampling_points_inner, sampling_points_outer);
else
    s_nm = compute_sample_sound_field_for_eq(c_nm, 0, 0, fs, taps_pw, N, grid_shape, normal_vector, c, rho, sampling_points);
end

% load MagLS-equalized HRTFs
hrtf_file_name = 'hrtfs/hrirs_ku100_%s_sh_N%d.mat';
hrirs_magls_sh = load(sprintf(hrtf_file_name, hrtf_type, N));

% -- get ground truth --
hrir_path = 'hrtfs/HRIR_L2702.sofa';
download_hrtfs(hrir_path);

SOFAstart;

hrirs_orig_sofa = SOFAload(hrir_path);

% --------------------- get EQ for left and right ear ---------------------
eq_l_ir = compute_dir_indep_sound_field_eq(s_nm, hrirs_magls_sh, hrirs_orig_sofa, fs, N, 1); 
eq_r_ir = compute_dir_indep_sound_field_eq(s_nm, hrirs_magls_sh, hrirs_orig_sofa, fs, N, 2); 

eq_ir = [eq_l_ir, eq_r_ir];

fprintf('done.\n\n');

% ----- verify auralization of anechoic data against the ground truth -----

fprintf('Computing anechoic auralization data for verification ... ');

fprintf('using head orientation %d deg ... ', round(head_orientation_azimuth_deg));

% s_nm is from a plane wave impinging from the front 

% test two different head orientations of the listener
brirs_0  = render_ambi_signals_binaurally_t(s_nm, hrirs_magls_sh, eq_ir,   (0+head_orientation_azimuth_deg)/180*pi, N); % sound incidence from straight ahead
brirs_90 = render_ambi_signals_binaurally_t(s_nm, hrirs_magls_sh, eq_ir, (-90+head_orientation_azimuth_deg)/180*pi, N); % sound incidence from the left

% this is used in the plots; it originates from the MagLS HRTFs computation
f_transition = max(1000, N*500); % Hz

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

output_file_name = sprintf('auralization_matrices/auralization_matrix_ambisonics_%s_%s_L%d.mat', data_type_string, grid_shape, size(sampling_points, 2));

fprintf('Storing the auralization matrix together with the sampling grid in file ''%s'' ... ', output_file_name);

if exist('sampling_points_outer', 'var')
    save(output_file_name, 'c_nm', 'fs', 'eq_ir', 'sampling_points_inner', 'sampling_points_outer', '-v7.3');
else
    save(output_file_name, 'c_nm', 'fs', 'eq_ir', 'sampling_points', '-v7.3');
end

fprintf('done.\n\n');
