function [pressure, velocity] = compute_plane_wave_on_grid(sampling_points, pw_incidence_azi_deg, pw_incidence_ele_deg, c, taps_pw, fs, normal_vector, rho)
%Simulates pressure and velocity due to a plane wave by means of a
% spherical harmonic synthesis
%
% The incidence direction is given by pw_incidence_azi_deg and 
% pw_incidence_ele_deg.
%
% normal_vector is only required if velocity is to be
% computed.
%
% taps_pw: 2-element vector (effective length of pw, length to zero-pad to)

if nargout == 2
    compute_velocity = 1;
    % assure unit length
    normal_vector = normal_vector ./ vecnorm(normal_vector);
else
    compute_velocity = 0;
    normal_vector = [];
end

% convert incidence direction to propagation direction
pw_prop_azi_deg =  pw_incidence_azi_deg + 180;
pw_prop_ele_deg = -pw_incidence_ele_deg;

% propagation direction in Cartesian coordinates
pw_prop_cart = [cos(pw_prop_azi_deg/180*pi) * sin(pi/2 - pw_prop_ele_deg/180*pi);
                sin(pw_prop_azi_deg/180*pi) * sin(pi/2 - pw_prop_ele_deg/180*pi);
                cos(pi/2 - pw_prop_ele_deg/180*pi)];

% Dot product to compute distance between sampling points and the wave
% front the moment the wave front passes the origin
% https://mathinsight.org/distance_point_plane
% https://se.mathworks.com/matlabcentral/answers/371665-distance-from-point-to-plane-plane-was-created-from-3d-point-data
distances = sampling_points.' * pw_prop_cart;

% convert to propagation time in samples
delay_samples = distances/c * fs;

% ----------------------- compute the pressure signal ---------------------
delay_samples_max = max(abs(delay_samples));
assert(taps_pw(1) > delay_samples_max(1));

all_taps = (0:taps_pw(1)/2).';

pressure = zeros(taps_pw(1)/2+1, size(sampling_points, 2));

if compute_velocity
    velocity = zeros(taps_pw(1)/2+1, size(sampling_points, 2));
end

% loop over all sampling points
for l = 1 : size(sampling_points, 2)
    
    % delay_samples is the offset to the middle of the buffer
    pressure(:, l) = exp(-1i .* all_taps/taps_pw(1) .* 2*pi .* (delay_samples(l)+taps_pw(1)/2));

    if compute_velocity
        velocity(:, l) = 1/(rho*c) .* pressure(:, l) .* dot(pw_prop_cart, normal_vector(:, l));
    end

end

pressure = ifft([pressure; conj(flipud(pressure(2:end-1, :)))], [], 1, 'symmetric');

if compute_velocity
    velocity = ifft([velocity; conj(flipud(velocity(2:end-1, :)))], [], 1, 'symmetric');
end

% ----------------------- window the irs, just in case --------------------
win      = hann(500); 
fade_in  = win(1:end/2);
fade_out = win(end/2+1:end);

pressure(1:length(fade_in),          :, :) = pressure(1:length(fade_in),          :, :) .* repmat(fade_in,  1, size(pressure, 2), size(pressure, 3));
pressure(end-length(fade_out)+1:end, :, :) = pressure(end-length(fade_out)+1:end, :, :) .* repmat(fade_out, 1, size(pressure, 2), size(pressure, 3));

if compute_velocity
    velocity(1:length(fade_in),          :, :) = velocity(1:length(fade_in),          :, :) .* repmat(fade_in,  1, size(velocity, 2), size(velocity, 3));
    velocity(end-length(fade_out)+1:end, :, :) = velocity(end-length(fade_out)+1:end, :, :) .* repmat(fade_out, 1, size(velocity, 2), size(velocity, 3));
end

% --------------------------------- zero pad ------------------------------
pressure = [pressure; zeros(taps_pw(2) - taps_pw(1), size(pressure, 2))];

if compute_velocity
    velocity = [velocity; zeros(taps_pw(2) - taps_pw(1), size(velocity, 2))];
end

end




