function [sampling_points_inner, sampling_points_outer] = create_grid_surface_cube(L, dx, delta)
% L: edge length
% dx: spacing of sampling points
% delta: distance between layers

spatial_axis = -L/2 + dx : dx : L/2 - dx;

% outer bottom layer
[x1, y1] = meshgrid(spatial_axis, spatial_axis);
z1 = -L/2 * ones(size(x1));

% inner bottom layer
x2 = x1;
y2 = y1;
z2 = z1 + delta;

sampling_points_outer    = [x1(:).'; y1(:).'; z1(:).'];
sampling_points_inner    = [x2(:).'; y2(:).'; z2(:).'];

% outer top layer
x3 = x1;
y3 = y1;
z3 = L/2 * ones(size(x3));

% inner top layer
x4 = x1;
y4 = y1;
z4 = z3 - delta;

sampling_points_outer    = [sampling_points_outer, [x3(:).'; y3(:).'; z3(:).']];
sampling_points_inner    = [sampling_points_inner, [x4(:).'; y4(:).'; z4(:).']];

% outer wall 1
[y5, z5] = meshgrid(spatial_axis, spatial_axis);
x5 = L/2 * ones(size(y5));

% inner wall 1
y6 = y5;
z6 = z5;
x6 = x5 - delta;

sampling_points_outer    = [sampling_points_outer, [x5(:).'; y5(:).'; z5(:).']];
sampling_points_inner    = [sampling_points_inner, [x6(:).'; y6(:).'; z6(:).']];

% outer wall 2
y7 = y5;
z7 = z5;
x7 = -L/2 * ones(size(y7));

% inner wall 2
y8 = y7;
z8 = z7;
x8 = x7 + delta;

sampling_points_outer    = [sampling_points_outer, [x7(:).'; y7(:).'; z7(:).']];
sampling_points_inner    = [sampling_points_inner, [x8(:).'; y8(:).'; z8(:).']];

% outer wall 3
[x9, z9] = meshgrid(spatial_axis, spatial_axis);
y9 = L/2 * ones(size(x9));

% inner wall 3
x10 = x9;
z10 = z9;
y10 = y9 - delta;

sampling_points_outer    = [sampling_points_outer, [x9(:).'; y9(:).'; z9(:).']];
sampling_points_inner    = [sampling_points_inner, [x10(:).'; y10(:).'; z10(:).']];

% outer wall 4
x11 = x9;
z11 = z9;
y11 = -L/2 * ones(size(x11));

% inner wall 4
x12 = x11;
z12 = z11;
y12 = y11 + delta;

sampling_points_outer = [sampling_points_outer, [x11(:).'; y11(:).'; z11(:).']];
sampling_points_inner = [sampling_points_inner, [x12(:).'; y12(:).'; z12(:).']];

% --------------------------- add edges -----------------------------------

% % top, parallel to x-axis
% x13 = [-L/2, spatial_axis, L/2];
% y13 = -L/2 * ones(size(x13));
% z13 =  L/2 * ones(size(x13));
% 
% x14 = x13;
% y14 = y13 + delta;
% z14 = z13 - delta;
% 
% grid.sampling_points_outer    = [grid.sampling_points_outer, [x13(:).'; y13(:).'; z13(:).']];
% grid.sampling_points_inner    = [grid.sampling_points_inner, [x14(:).'; y14(:).'; z14(:).']];
% 
% x15 = [-L/2, spatial_axis, L/2];
% y15 =  L/2 * ones(size(x15));
% z15 =  L/2 * ones(size(x15));
% 
% x16 = x15;
% y16 = y15 - delta;
% z16 = z15 - delta;
% 
% grid.sampling_points_outer    = [grid.sampling_points_outer, [x15(:).'; y15(:).'; z15(:).']];
% grid.sampling_points_inner    = [grid.sampling_points_inner, [x16(:).'; y16(:).'; z16(:).']];
% 
% % top, parallel to y-axis
% y17 = spatial_axis;
% x17 = -L/2 * ones(size(y17));
% z17 =  L/2 * ones(size(y17));
% 
% x18 = x17 + delta;
% y18 = y17;
% z18 = z17 - delta;
% 
% grid.sampling_points_outer    = [grid.sampling_points_outer, [x17(:).'; y17(:).'; z17(:).']];
% grid.sampling_points_inner    = [grid.sampling_points_inner, [x18(:).'; y18(:).'; z18(:).']];
% 
% y19 = spatial_axis;
% x19 =  L/2 * ones(size(y19));
% z19 =  L/2 * ones(size(y19));
% 
% x20 = x19 - delta;
% y20 = y19;
% z20 = z19 - delta;
% 
% grid.sampling_points_outer    = [grid.sampling_points_outer, [x19(:).'; y19(:).'; z19(:).']];
% grid.sampling_points_inner    = [grid.sampling_points_inner, [x20(:).'; y20(:).'; z20(:).']];



% figure;
% plot3(x1(:), y1(:), z1(:), '.');
% hold on;
% plot3(x2(:), y2(:), z2(:), '.');
% plot3(x3(:), y3(:), z3(:), '.');
% plot3(x4(:), y4(:), z4(:), '.');
% plot3(x5(:), y5(:), z5(:), '.');
% plot3(x6(:), y6(:), z6(:), '.');
% plot3(x7(:), y7(:), z7(:), '.');
% plot3(x8(:), y8(:), z8(:), '.');
% plot3(x9(:), y9(:), z9(:), '.');
% plot3(x10(:), y10(:), z10(:), '.');
% plot3(x11(:), y11(:), z11(:), '.');
% plot3(x12(:), y12(:), z12(:), '.');
% %plot3(x13(:), y13(:), z13(:), '.');
% hold off;
% axis equal;
% xlabel('x');
% ylabel('y');
% zlabel('z');

end

