function [x, y, z] = create_grid_surface_sphere(R, N)
%Creates a cubical volumetric sampling grid
%
% D:         Side length of the cube in m
% L:         Number of sampling points per dimension
% N_desired: Desired SH order to be computed _from_ the data; is used for 
%            predicting required no. of points
%
% The sound field is always simulated at 32nd order.
% The incidence elevation is always 0 deg.

% ------ create spherical mesh of sampling points with two layers ---------

% --- Lebedev ---
% % %grid_tmp = load('resources/lebedev_grid_2702.mat');
% grid_tmp = load('resources/lebedev_grids/lebedev_grid_110.mat');
% grid_tmp.ele = pi/2 - grid_tmp.col;
% % 
% [x, y, z] = sph2cart(grid_tmp.azi, grid_tmp.ele, R); 

%  --- Fliege ---
grid_tmp = load(sprintf('resources/fliege_grids/%d.txt', (N+1)^2));

x = R * grid_tmp(:, 1).';
y = R * grid_tmp(:, 2).';
z = R * grid_tmp(:, 3).';


% --- t-design ---
%grid_tmp = load('resources/t_designs/t_design_deg_10.mat'); 
%grid_tmp = load('resources/t_designs/t_design_deg_16.mat'); % high resolution

% x = R * grid_tmp.vecs(:, 1).';
% y = R * grid_tmp.vecs(:, 2).';
% z = R * grid_tmp.vecs(:, 3).';




