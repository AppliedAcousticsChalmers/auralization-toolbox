function [sampling_points_inner, sampling_points_outer] = get_grid_surface_cube(L, dx, delta)
% L: edge length
% dx: spacing of sampling points
% delta: distance between layers
%
% The outer layer has an edge length of precisely L. The edge length of the
% inner layer is L-2*delta.

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

% top, parallel to x-axis
x13 = [-L/2, spatial_axis, L/2];
y13 = -L/2 * ones(size(x13));
z13 =  L/2 * ones(size(x13));

x14      = x13;
x14(1)   = x14(1)   + delta;
x14(end) = x14(end) - delta;
y14 = y13 + delta;
z14 = z13 - delta;

sampling_points_outer    = [sampling_points_outer, [x13(:).'; y13(:).'; z13(:).']];
sampling_points_inner    = [sampling_points_inner, [x14(:).'; y14(:).'; z14(:).']];

x15 = [-L/2, spatial_axis, L/2];
y15 =  L/2 * ones(size(x15));
z15 =  L/2 * ones(size(x15));

x16      = x15;
x16(1)   = x16(1)   + delta;
x16(end) = x16(end) - delta;
y16 = y15 - delta;
z16 = z15 - delta;

sampling_points_outer    = [sampling_points_outer, [x15(:).'; y15(:).'; z15(:).']];
sampling_points_inner    = [sampling_points_inner, [x16(:).'; y16(:).'; z16(:).']];

% top, parallel to y-axis
y17 = spatial_axis;
x17 = -L/2 * ones(size(y17));
z17 =  L/2 * ones(size(y17));

x18 = x17 + delta;
y18 = y17;
z18 = z17 - delta;

sampling_points_outer    = [sampling_points_outer, [x17(:).'; y17(:).'; z17(:).']];
sampling_points_inner    = [sampling_points_inner, [x18(:).'; y18(:).'; z18(:).']];

y19 = spatial_axis;
x19 =  L/2 * ones(size(y19));
z19 =  L/2 * ones(size(y19));

x20 = x19 - delta;
y20 = y19;
z20 = z19 - delta;

sampling_points_outer    = [sampling_points_outer, [x19(:).'; y19(:).'; z19(:).']];
sampling_points_inner    = [sampling_points_inner, [x20(:).'; y20(:).'; z20(:).']];

% bottom, parallel to x-axis
x21 = [-L/2, spatial_axis, L/2];
y21 = -L/2 * ones(size(x21));
z21 = -L/2 * ones(size(x21));

x22      = x21;
x22(1)   = x22(1)   + delta;
x22(end) = x22(end) - delta;
y22 = y21 + delta;
z22 = z21 + delta;

sampling_points_outer    = [sampling_points_outer, [x21(:).'; y21(:).'; z21(:).']];
sampling_points_inner    = [sampling_points_inner, [x22(:).'; y22(:).'; z22(:).']];

x23 =  x21;
y23 = -y21;
z23 =  z21;

x24      = x23;
x24(1)   = x24(1)   + delta;
x24(end) = x24(end) - delta;
y24 = y23 - delta;
z24 = z23 + delta;

sampling_points_outer    = [sampling_points_outer, [x23(:).'; y23(:).'; z23(:).']];
sampling_points_inner    = [sampling_points_inner, [x24(:).'; y24(:).'; z24(:).']];

% bottom, parallel to y-axis
y25 = spatial_axis;
x25 = -L/2 * ones(size(y17));
z25 = -L/2 * ones(size(y17));

x26 = x25 + delta;
y26 = y25;
z26 = z25 + delta;

sampling_points_outer    = [sampling_points_outer, [x25(:).'; y25(:).'; z25(:).']];
sampling_points_inner    = [sampling_points_inner, [x26(:).'; y26(:).'; z26(:).']];

y27 = spatial_axis;
x27 =  L/2 * ones(size(y19));
z27 = -L/2 * ones(size(y19));

x28 = x27 - delta;
y28 = y27;
z28 = z27 + delta;

sampling_points_outer    = [sampling_points_outer, [x27(:).'; y27(:).'; z27(:).']];
sampling_points_inner    = [sampling_points_inner, [x28(:).'; y28(:).'; z28(:).']];

% add vertical edge (x, y) = (L/2, -L/2)
x29 =  L/2 * ones(size(spatial_axis)); 
y29 = -L/2 * ones(size(spatial_axis));
z29 =  spatial_axis;

x30 = x29 - delta;
y30 = y29 + delta;
z30 = z29;

sampling_points_outer    = [sampling_points_outer, [x29(:).'; y29(:).'; z29(:).']];
sampling_points_inner    = [sampling_points_inner, [x30(:).'; y30(:).'; z30(:).']];

% add vertical edge (x, y) = (L/2, L/2)
x31 =  L/2 * ones(size(spatial_axis)); 
y31 =  L/2 * ones(size(spatial_axis));
z31 =  spatial_axis;

x32 = x31 - delta;
y32 = y31 - delta;
z32 = z31;

sampling_points_outer    = [sampling_points_outer, [x31(:).'; y31(:).'; z31(:).']];
sampling_points_inner    = [sampling_points_inner, [x32(:).'; y32(:).'; z32(:).']];

% add vertical edge (x, y) = (-L/2, L/2)
x33 = -L/2 * ones(size(spatial_axis)); 
y33 =  L/2 * ones(size(spatial_axis));
z33 =  spatial_axis;

x34 = x33 + delta;
y34 = y33 - delta;
z34 = z33;

sampling_points_outer    = [sampling_points_outer, [x33(:).'; y33(:).'; z33(:).']];
sampling_points_inner    = [sampling_points_inner, [x34(:).'; y34(:).'; z34(:).']];

% add vertical edge (x, y) = (-L/2, -L/2)
x35 = -L/2 * ones(size(spatial_axis)); 
y35 = -L/2 * ones(size(spatial_axis));
z35 =  spatial_axis;

x36 = x35 + delta;
y36 = y35 + delta;
z36 = z35;

sampling_points_outer    = [sampling_points_outer, [x35(:).'; y35(:).'; z35(:).']];
sampling_points_inner    = [sampling_points_inner, [x36(:).'; y36(:).'; z36(:).']];

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
% plot3(x13(:), y13(:), z13(:), '.');
% plot3(x14(:), y14(:), z14(:), '.');
% plot3(x15(:), y15(:), z15(:), '.');
% plot3(x16(:), y16(:), z16(:), '.');
% plot3(x17(:), y17(:), z17(:), '.');
% plot3(x18(:), y18(:), z18(:), '.');
% plot3(x19(:), y19(:), z19(:), '.');
% plot3(x20(:), y20(:), z20(:), '.');
% plot3(x21(:), y21(:), z21(:), '.');
% plot3(x22(:), y22(:), z22(:), '.');
% plot3(x23(:), y23(:), z23(:), '.');
% plot3(x24(:), y24(:), z24(:), '.');
% plot3(x25(:), y25(:), z25(:), '.');
% plot3(x26(:), y26(:), z26(:), '.');
% plot3(x27(:), y27(:), z27(:), '.');
% plot3(x28(:), y28(:), z28(:), '.');
% plot3(x29(:), y29(:), z29(:), '.');
% plot3(x30(:), y30(:), z30(:), '.');
% plot3(x31(:), y31(:), z31(:), '.');
% plot3(x32(:), y32(:), z32(:), '.');
% plot3(x33(:), y33(:), z33(:), '.');
% plot3(x34(:), y34(:), z34(:), '.');
% plot3(x35(:), y35(:), z35(:), '.');
% plot3(x36(:), y36(:), z36(:), '.');
% hold off;
% axis equal;
% xlabel('x');
% ylabel('y');
% zlabel('z');

end

