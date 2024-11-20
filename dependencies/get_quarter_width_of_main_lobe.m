function [quarter_width_rad] = get_quarter_width_of_main_lobe(N)

% empirical data
half_width_rad = [180 109 72.7 55.5 43.6 32 24.5 21.8 20 18 12 10.5 7  2.17]/180*pi;
order          = [ 0   1    2    3   4    6  8    9   10 11 17  20  30  100];

% get quarter width of main lobe for the given order
quarter_width_rad = interp1(order, half_width_rad/2, N); 



end

