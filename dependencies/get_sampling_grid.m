function [output_1, output_2, normal_vector] = get_sampling_grid(grid_shape, layer_type, R, L_cube, L_sphere, delta)
% Output:
%    output_1: sampling points (if single layer), 
%              inner-layer sampling points (if double layer)
%    output_2: surface normal vector (if single layer), 
%              outer-layer sampling points (if double layer)
%
% Input: 
%     grid_shape: 'cubical_volume', 'spherical_surface', or 'cubical_surface'
%     layer_type: 'single' or 'double'
%     R:          Radius (if spherical) or half edge length (if cubical) of grid
%     L_cube:     Number of sampling points per dimension (if cubical grid)
%     L_sphere:   Total number of sampling points (if spherical grid);
%                                must be a square number in this case
%     delta:      Distance between layers in a double layer grid
%
% Example for a double layer: 
% [sampling_points_inner, sampling_points_outer, normal_vector] = get_sampling_grid('spherical_surface', 'double', R, L, delta);
% Example for a single layer: 
% [sampling_points,              ~             , normal_vector] = get_sampling_grid('spherical_surface', 'single', R, L);
%

if strcmp(grid_shape, 'cubical_volume')
    
    sampling_points = get_grid_volumetric_cube(2*R, L_cube);
    
    % make it a ball
    %r = sqrt(x.^2 + y.^2 + z.^2);
    %sampling_points = sampling_points(:, r <= D*.6);
    
    %plot_grid();

    output_1      = sampling_points;
    output_2      = [];
    normal_vector = [];

elseif strcmp(grid_shape, 'spherical_surface')
        
    if strcmp(layer_type, 'single') 

        sampling_points = get_grid_surface_sphere(R, L_sphere);

        % surface normal
        normal_vector = sampling_points/R;
       
        output_1      = sampling_points;
        output_2      = [];

    elseif strcmp(layer_type, 'double')

        sampling_points_outer = get_grid_surface_sphere(R+delta/2, L_sphere);
        sampling_points_inner = get_grid_surface_sphere(R-delta/2, L_sphere);
  
        % surface normal
        normal_vector = sampling_points_outer - sampling_points_inner;
        normal_vector = normal_vector ./ vecnorm(normal_vector);

        % for file storage and to simplify the syntax
        %sampling_points = (sampling_points_inner + sampling_points_outer)/2; 
        
        output_1      = sampling_points_inner;
        output_2      = sampling_points_outer;

    end
    
elseif strcmp(grid_shape, 'cubical_surface')
                
    [sampling_points_inner, sampling_points_outer] = get_grid_surface_cube(2*R, 2*R/(L_cube-1), delta);
    
    % surface normal
    normal_vector = sampling_points_outer - sampling_points_inner;
    normal_vector = normal_vector ./ vecnorm(normal_vector);
   
    % cubical single-layer surface
    if strcmp(layer_type, 'single') 

        output_1 = sampling_points_outer;
        output_2 = [];
        
    else
        output_1 = sampling_points_inner;
        output_2 = sampling_points_outer;
    end
    
else
    error('Unknown grid shape.');  

end % if cubical or spherical sampling

end
