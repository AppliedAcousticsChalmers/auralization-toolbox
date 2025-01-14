function [s_nm, sampled_sound_field] = compute_sample_sound_field_for_eq(c_nm, azi_pw_deg, ele_pw_deg, fs, taps_pw, N, grid_shape, normal_vector, c, rho, varargin)
% Computes the ambisonic signal for a plane wave that impinges from
% straight head onto the grid. The sound field equalization in computed
% based on this.

if length(varargin) == 1
    
    sampling_points = varargin{1};
    layer_type = 'single'; 

elseif length(varargin) == 2

    sampling_points_inner = varargin{1}; 
    sampling_points_outer = varargin{2};
    layer_type = 'double'; 

else
    error('Too many input arguments.');

end

% ---------------------- compute sound field on grid ----------------------
% ----------------------    plane wave from (0, 0)   ----------------------

if strcmp(grid_shape, 'cubical_volume')
        
    pressure = compute_plane_wave_on_grid(sampling_points, azi_pw_deg, ele_pw_deg, c, taps_pw, fs);
    
elseif strcmp(grid_shape, 'spherical_surface')
        
    if strcmp(layer_type, 'single') 
        [pressure, velocity] = compute_plane_wave_on_grid(sampling_points, azi_pw_deg, ele_pw_deg, c, taps_pw, fs, normal_vector, rho);
    elseif strcmp(layer_type, 'double')
        pressure_outer = compute_plane_wave_on_grid(sampling_points_outer, azi_pw_deg, ele_pw_deg, c, taps_pw, fs, normal_vector, rho);
        pressure_inner = compute_plane_wave_on_grid(sampling_points_inner, azi_pw_deg, ele_pw_deg, c, taps_pw, fs, normal_vector, rho);
    end

elseif strcmp(grid_shape, 'cubical_surface')
     
    % cubical single-layer surface
    if strcmp(layer_type, 'single') 

        [pressure, velocity] = compute_plane_wave_on_grid(sampling_points, azi_pw_deg, ele_pw_deg, c, taps_pw, fs, normal_vector, rho);

    elseif strcmp(layer_type, 'double')
    
        pressure_outer = compute_plane_wave_on_grid(sampling_points_outer, azi_pw_deg, ele_pw_deg, c, taps_pw, fs, normal_vector, rho);
        pressure_inner = compute_plane_wave_on_grid(sampling_points_inner, azi_pw_deg, ele_pw_deg, c, taps_pw, fs, normal_vector, rho);

    end

end

% --------------------- compute virtual cardioids -------------------------

% prepare the sound field data for SH decomposition
if strcmp(grid_shape, 'spherical_surface') || strcmp(grid_shape, 'cubical_surface')
    
    if exist('velocity', 'var')
        % compute the signals captured by the virtual cardioid transducers
        sampled_sound_field = pressure - rho*c .* velocity;
    else
        % compute the signals captured by the virtual cardioid transducers
        [~, sampled_sound_field] = compute_cardioid_from_pressure(sampling_points_inner, sampling_points_outer, pressure_inner, pressure_outer, fs, c);
    end
    
else % if volumetric
    
    sampled_sound_field = pressure;
    
end

% window
%win = hann(100);
%fade_in = win(1:end/2);
%fade_out = win(end/2+1:end);

%sampled_sound_field(         1:length(fade_in), :) = sampled_sound_field(1:length(fade_in)         , :) .* repmat(fade_in,  1, size(sampled_sound_field, 2));
%sampled_sound_field(end-length(fade_out)+1:end, :) = sampled_sound_field(end-length(fade_out)+1:end, :) .* repmat(fade_out, 1, size(sampled_sound_field, 2));


% ------ if ambisonic representation of the room data is desired ----------
if size(c_nm, 1) > 1
    s_nm = get_sh_coefficients_t(sampled_sound_field, c_nm, N);
else
    s_nm = [];
end


% zero pad
%s_nm = [s_nm; zeros(2048, size(s_nm, 2))];

% figure;
% plot(20*log10(abs(sampled_sound_field)))
% 
% figure;
% plot(20*log10(abs(c_nm(:, :, 2))))
% 
% 
% figure;
% plot(20*log10(abs(s_nm) + 0.0001))

end

