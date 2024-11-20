function [out_binaural] = render_ambi_signals_binaurally_t(ambi_signals, hrirs, eq_ir, head_orientation_azimuth_rad, N)
% So far, this function works only for sphharm_type = 'real'. It is not
% conventient to extend it for complex SHs because this would produce
% time-domain output signals that are complex.
%
% The equation numbers refer to 
% 
%   Jens Ahrens, "Ambisonic Encoding of Signals From Spherical Microphone
%   Arrays," Technical note v. 1, Chalmers University of Technology, 2022.
%   https://arxiv.org/abs/2211.00583
%
% Written by Jens Ahrens, 2022

% ------------------------- Rendering, Eq. (12) ---------------------------

% zero pad for the convolution
ambi_signals = [ambi_signals; zeros(size(hrirs.h_l_ring, 1), size(ambi_signals, 2))];

out_binaural = zeros(size(ambi_signals, 1), 2);

for n = 0 : N
    for m = -n : n

        out_binaural(:, 1) = out_binaural(:, 1) + fftfilt(hrirs.h_l_ring(:, n^2+n+m+1) .* cos(m*head_orientation_azimuth_rad) - hrirs.h_l_ring(:, n^2+n-m+1) .* sin(m*head_orientation_azimuth_rad), ambi_signals(:, n^2+n+m+1));
        out_binaural(:, 2) = out_binaural(:, 2) + fftfilt(hrirs.h_r_ring(:, n^2+n+m+1) .* cos(m*head_orientation_azimuth_rad) - hrirs.h_r_ring(:, n^2+n-m+1) .* sin(m*head_orientation_azimuth_rad), ambi_signals(:, n^2+n+m+1));

    end
end

% apply EQ
if ~isempty(eq_ir)
    out_binaural = [out_binaural; zeros(size(eq_ir, 1), 2)];
    out_binaural = fftfilt(eq_ir, out_binaural);
end

end
