function [hrirs_rot_sofa] = rotate_ku100_hrirs(hrirs_sofa, azimuth_deg)
% rotates the acoustic data as if the head was oriented in azimuth_deg
% during the measurement (instead in 0°)

azimuth_rad = azimuth_deg/180*pi;

% extract irs
hrirs_l = squeeze(hrirs_sofa.Data.IR(:, 1, :)).';
hrirs_r = squeeze(hrirs_sofa.Data.IR(:, 2, :)).';

% extract incidence angles
azi_hrtfs =        hrirs_sofa.SourcePosition(:, 1).' / 180 * pi;
col_hrtfs = pi/2 - hrirs_sofa.SourcePosition(:, 2).' / 180 * pi; % ele to col

% Get Lebedev grid weights for the quadrature of the transformation 
% integral; not need if least squares SH fit is performed.
tmp = load('resources/2702_lebedev_grid.mat');
grid_weights = tmp.grid_weights;
N = tmp.max_order;
clear tmp;

h_l_ring = zeros(size(hrirs_l, 1), (N+1)^2);
h_r_ring = zeros(size(hrirs_r, 1), (N+1)^2);

% --- SH decomposition ---
for n = 0 : N 
    for m = -n : n        
        h_l_ring(:, n^2+n+m+1) = sum(hrirs_l .* repmat(grid_weights .* sphharm(n, m, col_hrtfs, azi_hrtfs, 'real'), [size(hrirs_l, 1) 1]), 2) * 4*pi;
        h_r_ring(:, n^2+n+m+1) = sum(hrirs_r .* repmat(grid_weights .* sphharm(n, m, col_hrtfs, azi_hrtfs, 'real'), [size(hrirs_r, 1) 1]), 2) * 4*pi;
    end
end

hrirs_l_rot = zeros(size(hrirs_l));
hrirs_r_rot = zeros(size(hrirs_r));

% --- SH recomposition ---
for n = 0 : N 
    for m = -n : n        
        hrirs_l_rot = hrirs_l_rot + h_l_ring(:, n^2+n+m+1) .* sphharm(n, m, col_hrtfs, azi_hrtfs - azimuth_rad, 'real');
        hrirs_r_rot = hrirs_r_rot + h_r_ring(:, n^2+n+m+1) .* sphharm(n, m, col_hrtfs, azi_hrtfs - azimuth_rad, 'real');
    end
end

% return it in SOFA format
hrirs_rot              = cat(3, hrirs_l_rot, hrirs_r_rot);
hrirs_rot_sofa         = hrirs_sofa;
hrirs_rot_sofa.Data.IR = permute(hrirs_rot, [2, 3, 1]);

%hrirs_rot_sofa.SourcePosition(:, 1) = hrirs_sofa.SourcePosition(:, 1) - azimuth_deg;

end