clear;

addpath('dependencies/');

%grid_shape = 'cubical_volume'; L_cube = 5;
grid_shape = 'spherical_surface'; L_sphere = 81;
%grid_shape = 'cubical_surface'; L_cube = 5;

% L_cube is the number of sampling points per dimension (if cubical grid)
% L_sphere is the total number of sampling points (if spherical grid);    
%                               must be a square number in this case

layer_type = 'single'; % pressure and velocity are computed for surface grids
%layer_type = 'double'; % only pressure is computed

R     = .07;  % radius of sphere / half edge length of cube
delta = .001; % distance between layers (if double layer, disregrded otherwise)

output_file_name = 'resources/grid.mat';

% ------------------------------------------------------------------------
% L_cube | max. N | no. of sampling points for cubical volumetric grids:
% ------------------------
%    3   |    2   |   27
%    4   |    3   |   64
%    5   |    5   |  125
%    6   |    7   |  216 
%    7   |   11   |  343
%    8   |   15   |  512
%    9   |   18   |  729
%   10   |   18   | 1000
%   11   |   18   | 1331
%   12   |   18   | 1728
%   13   |   18   | 2197
%
% L_cube | max. N | no. of sampling points for cubical surface grids:
% ------------------------
%    2   |    3   |    26
%    3   |    5   |    56
%    4   |    7   |    98
%    5   |    9   |   152
%    6   |   11   |   218
%    7   |   13   |   296
%    8   |   15   |   386
%    9   |   17   |   488
%   10   |   19   |   602
% ------------------------------------------------------------------------

% ------------------ avoid confusion and syntax errors -------------------

if strcmp(grid_shape, 'cubical_volume')
    layer_type = '';
    L_sphere = NaN;
end

if strcmp(grid_shape, 'spherical_surface')
    L_cube = NaN;
end

if strcmp(grid_shape, 'cubical_surface')
    L_sphere = NaN;
end

% ---------------------------- get grid ----------------------------------


% ---- get the spatial grid on which to compute the pressure/velocity ----
[output_1, output_2, normal_vector] = get_sampling_grid(grid_shape, layer_type, R, L_cube, L_sphere, delta);

% sanity check
if isnan(output_1)
    error('Something''s wrong here.');
end

% sort the output data
if ~strcmp(grid_shape, 'cubical_volume') && strcmp(layer_type, 'double')
    sampling_points_inner = output_1;
    sampling_points_outer = output_2;

    % combine inner and outer layer in one variable for simpler implementation
    sampling_points = [sampling_points_inner, sampling_points_outer];
else
    sampling_points = output_1;
    % output_2 is not used
end

% ---------------------------- store grid ---------------------------------

if strcmp(grid_shape, 'cubical_volume') || strcmp(layer_type, 'single')
    save(output_file_name, 'sampling_points', 'grid_shape', 'layer_type', '-v7.3');
elseif strcmp(layer_type, 'double')
    save(output_file_name, 'sampling_points_inner', 'sampling_points_outer', 'grid_shape', 'layer_type', '-v7.3');
else
    error('Unknown combination of grid_shape an layer_type.');
end

fprintf('\nStored ''%s'' grid in file ''%s''.\n\n', grid_shape, output_file_name);

% --------------------------- plot grid -----------------------------------

plot_sampling_grid(output_file_name);

