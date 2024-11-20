function [sampling_points] = get_grid_surface_sphere(R, L)
%Loads a Fliege grid with L points
%
% R:         Radius of the sphere in m

% ------ create spherical mesh of sampling points with two layers ---------

% --- Lebedev ---
% % %grid_tmp = load('resources/lebedev_grid_2702.mat');
% grid_tmp = load('resources/lebedev_grids/lebedev_grid_110.mat');
% grid_tmp.ele = pi/2 - grid_tmp.col;
% % 
% [x, y, z] = sph2cart(grid_tmp.azi, grid_tmp.ele, R); 

%  --- Fliege ---
grid_tmp = load(sprintf('resources/fliege_grids/%d.txt', L));

%addpath('/Users/jensah/Documents/coding/polarch/Spherical-Harmonic-Transform/');
%grid_tmp = getTdesign(N);

x = R * grid_tmp(:, 1).';
y = R * grid_tmp(:, 2).';
z = R * grid_tmp(:, 3).';

% --- t-design ---
%grid_tmp = load('resources/t_designs/t_design_deg_10.mat'); 
%grid_tmp = load('resources/t_designs/t_design_deg_16.mat'); % high resolution

% x = R * grid_tmp.vecs(:, 1).';
% y = R * grid_tmp.vecs(:, 2).';
% z = R * grid_tmp.vecs(:, 3).';

sampling_points = [x; y; z];

end
