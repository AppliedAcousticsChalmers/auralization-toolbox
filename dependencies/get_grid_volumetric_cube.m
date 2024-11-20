function [sampling_points] = get_grid_volumetric_cube(D, L)
%Creates a cubical volumetric sampling grid
%
% D:         Side length of the cube in m
% L:         Number of sampling points per dimension
% N_desired: Desired SH order to be computed _from_ the data; is used for 
%            predicting required no. of points
%
% The sound field is always simulated at 32nd order.
% The incidence elevation is always 0 deg.

dx = D/(L-1); % sampling interval
 
% grid_data.grid_type = 'volume';

% --------- create mesh of sampling points; cube shaped volume -----------

spatial_axis = -D/2 : dx : D/2; % m
[x, y, z] = meshgrid(spatial_axis, spatial_axis, spatial_axis);

% reshape
x = x(:).';
y = y(:).';
z = z(:).';

sampling_points = [x; y; z];



