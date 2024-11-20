function plot_sampling_grid(vargin)

% Either call the script as a function or change the following default value of 
% 'grid_file' to run it as a self-suficient script.

if nargin ~= 1
    % sampling grid
    grid_file = 'resources/grid_spherical_surface_L25.mat';
else
    grid_file = vargin;
end

% -------------------------------------------------------------------------

load(grid_file);

figure;
set(gcf, 'Color', [1 1 1]);
set(gcf, 'Position', [700 100 500 500]);

if strcmp(layer_type, 'double')

    plot3(sampling_points_outer(1, :), sampling_points_outer(2, :), sampling_points_outer(3, :), '.');
    hold on;
    plot3(sampling_points_inner(1, :), sampling_points_inner(2, :), sampling_points_inner(3, :), '.');
    hold off;

else
    plot3(sampling_points(1, :), sampling_points(2, :), sampling_points(3, :), '.', 'MarkerSize', 12);
end

grid on;
box on;
axis equal;

% find farthest sampling point
r = vecnorm(sampling_points, 2, 1);
axis([-1 1 -1 1 -1 1] * max(r)*1.2);

xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');

camproj('perspective');

end
