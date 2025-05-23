function [c_l, c_r, condition_number] = compute_c_direct(head_orientation_deg, hrirs_sofa, azi_fliege_deg, ele_fliege_deg, c, taps_c_nm, taps_pw, f_transition, fs, dynamic_range_dB, grid_shape, layer_type, normal_vector, rho, sampling_points, varargin)
% rotates the HRTF set to account for the desired head orientation
% computes anechoic data for the least-squares fit and performs that fit

if length(varargin) == 2
    sampling_points_inner = varargin{1}; 
    sampling_points_outer = varargin{2};
end

% --------------- compute sample plane waves for the LS fit ---------------

fprintf('Rotating HRTF set to account for head orientation ... ');

% rotate HRIRs to account for the head orientation
hrirs_rot_sofa = rotate_ku100_hrirs(hrirs_sofa, head_orientation_deg);

fprintf('done.\n\n');

% pick the HRTFs for the least-squares fit
[azi_pw_deg, ele_pw_deg, indices_hrirs_calibration] = pick_angles_according_to_hrtfs(hrirs_rot_sofa, azi_fliege_deg, ele_fliege_deg);

% dimensions:                 taps  x  no. of sampling points  x  no. of plane waves
sampled_sound_field_tmp = zeros(taps_c_nm, size(sampling_points, 2), length(azi_pw_deg));

% -------------------- compute the calibration data -----------------------
display_progress('Computing the calibration data:');

for index = 1 : length(azi_pw_deg)

    display_progress(index/length(azi_pw_deg));

    if strcmp(grid_shape, 'cubical_volume') 

        pressure_tmp = compute_plane_wave_on_grid(sampling_points, azi_pw_deg(index), ele_pw_deg(index), c, [taps_pw(1) taps_c_nm], fs, normal_vector, rho);

        % compute the signals captured by the virtual cardioid transducers
        sampled_sound_field_tmp(:, :, index) = pressure_tmp;

    %if exist('velocity', 'var')
    elseif strcmp(layer_type, 'single')

        [pressure_tmp, velocity_tmp] = compute_plane_wave_on_grid(sampling_points, azi_pw_deg(index), ele_pw_deg(index), c, [taps_pw(1) taps_c_nm], fs, normal_vector, rho);

        % compute the signals captured by the virtual cardioid transducers
        sampled_sound_field_tmp(:, :, index) = pressure_tmp - rho*c .* velocity_tmp;

    elseif strcmp(layer_type, 'double')

        pressure_outer_tmp = compute_plane_wave_on_grid(sampling_points_outer, azi_pw_deg(index), ele_pw_deg(index), c, [taps_pw(1) taps_c_nm], fs);
        pressure_inner_tmp = compute_plane_wave_on_grid(sampling_points_inner, azi_pw_deg(index), ele_pw_deg(index), c, [taps_pw(1) taps_c_nm], fs);

        % TODO: Implement computation of particle velocity
        velocity_tmp = compute_velocity_from_double_pressure(sampling_points_inner, sampling_points_outer, pressure_inner_tmp, pressure_outer_tmp, fs, rho);
        
        sound_field_tmp = pressure_outer_tmp - rho*c .* velocity_tmp;

        sampled_sound_field_tmp(:, :, index) = sound_field_tmp;

    end

end % loop over plane wave incidence directions

fprintf('\n\n');

% --------------------- compute quadrature matrix -------------------------
     
[c_l, c_r, condition_number] = get_c_direct(sampled_sound_field_tmp, squeeze(hrirs_rot_sofa.Data.IR(indices_hrirs_calibration, 1, :)).', squeeze(hrirs_rot_sofa.Data.IR(indices_hrirs_calibration, 2, :)).', fs, f_transition, dynamic_range_dB);

end
