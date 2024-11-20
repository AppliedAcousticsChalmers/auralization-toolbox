function [points] = points_on_sphere(N)
% after Markus Deserno, "How to generate equidistributed points of the
% surface of a sphere"
% radius of sphere is 1

r = 1;

a = 4*pi*r^2/N;
d = sqrt(a);

M_col = round(pi/d);
d_col = pi/M_col; 
d_azi = a/d_col;

points = zeros(3, 0);

% loop over colatitude
for m = 0 : M_col-1

    col = pi*(m+0.5)/M_col;

    M_azi = round(2*pi*sin(col)/d_col);

    % loop over azimuth
    for n = 0 : M_azi-1

        azi    = 2*pi*n/M_azi;
        point  = [r*cos(azi)*sin(col); r*sin(azi)*sin(col); r*cos(col)];
        points = [points, point];

    end

end

% figure;
% plot3(points(1, :), points(2, :), points(3, :), '.', 'Markersize', 20);
% grid on;
% axis square;

end

