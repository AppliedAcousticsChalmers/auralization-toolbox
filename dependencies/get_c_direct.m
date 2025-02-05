function [c_l_irs, c_r_irs, condition_number] = get_c_direct(irs_calibration, hrirs_calibration_l, hrirs_calibration_r, fs, f_transition, dynamic_range_dB)
%
% irs_calibration: input signal
% hrirs_calibration_l, hrirs_calibration_r: desired output signal

%f_transition = 40000;

f = linspace(0, fs/2, size(irs_calibration, 1)/2+1).';

grp_delay_irs_s = zeros(0, 1);

% --- remove group delay from calibration data --- 
for ch = 1 : size(irs_calibration, 2)
    grp_delay_irs_s = [grp_delay_irs_s; grpdelay(irs_calibration(:, ch, 1), 1, f, fs)]; % only do it for one incidence direction
end

grp_delay_irs = median(grp_delay_irs_s);

% remove group delay
for dir = 1 : size(irs_calibration, 3)
    irs_calibration(:, :, dir) = fractional_delay(irs_calibration(:, :, dir), -grp_delay_irs);
end

% DFT
tfs = fft(irs_calibration, [], 1);
tfs = tfs(1:end/2+1, :, :);

% zero pad HRIRs
hrirs_calibration_l = [hrirs_calibration_l; ...
    zeros(size(irs_calibration, 1) - size(hrirs_calibration_l, 1), size(hrirs_calibration_l, 2))];
hrirs_calibration_r = [hrirs_calibration_r; ...
    zeros(size(irs_calibration, 1) - size(hrirs_calibration_r, 1), size(hrirs_calibration_r, 2))];

grp_delay_hrirs_l_s = zeros(0, 1);
grp_delay_hrirs_r_s = zeros(0, 1);

% --- remove group delay from HRTFs ---
for ch = 1 : size(hrirs_calibration_l, 2)
    grp_delay_hrirs_l_s = [grp_delay_hrirs_l_s; grpdelay(hrirs_calibration_l(:, ch), 1, f, fs)];
    grp_delay_hrirs_r_s = [grp_delay_hrirs_r_s; grpdelay(hrirs_calibration_r(:, ch), 1, f, fs)];
end

grp_delay_hrirs_l   = median(grp_delay_hrirs_l_s);
grp_delay_hrirs_r   = median(grp_delay_hrirs_r_s);
hrirs_calibration_l = fractional_delay(hrirs_calibration_l, -grp_delay_hrirs_l);
hrirs_calibration_r = fractional_delay(hrirs_calibration_r, -grp_delay_hrirs_r);

hrtfs_l = fft(hrirs_calibration_l, [], 1);
hrtfs_l = hrtfs_l(1:end/2+1, :);

hrtfs_r = fft(hrirs_calibration_r, [], 1);
hrtfs_r = hrtfs_r(1:end/2+1, :);

condition_number = zeros(size(tfs, 1), 1);

C_l = zeros(size(tfs, 1), size(tfs, 2));
C_r = zeros(size(tfs, 1), size(tfs, 2));

% -------------------------------------------------------------------------

% loop over all bins
for bin = 1 : size(tfs, 1)

    if f(bin) < 100 
        % this helps much with cardioid nodes
        svd_reg = 10^(-dynamic_range_dB(1)/20); % 20 dB
        %svd_reg = 0.03; % 30 dB (should be fine, too; to be confirmed)
    else
        svd_reg = 10^(-dynamic_range_dB(2)/20); % 40 dB
    end

    pw_grid = squeeze(tfs(bin, :, :));

    %[U, s, V] = svd(pw_grid.', 'econ', 'vector');
    [U, s, V] = svd(pw_grid.', 'econ'); s = diag(s); % for older MATLAB versions
    s = 1 ./ max(s, svd_reg * max(s)); % regularize
    Y_reg_inv = conj(U) * (s .* V.');

    if f(bin) < f_transition % least-squares below 

        C_l(bin, :) = hrtfs_l(bin, :) * Y_reg_inv;
        C_r(bin, :) = hrtfs_r(bin, :) * Y_reg_inv;

    elseif f(bin) < fs/2*0.85 % magnitude least-squares above 

        phi_l = angle(C_l(bin-1, :) * pw_grid);
        phi_r = angle(C_r(bin-1, :) * pw_grid);

        if bin == size(tfs, 1) && ~mod(size(irs_calibration, 1), 2) % Nyquist bin, is even
           C_l(bin, :) = real(abs(hrtfs_l(bin, :)) .* exp(1i * phi_l)) * Y_reg_inv;
           C_r(bin, :) = real(abs(hrtfs_r(bin, :)) .* exp(1i * phi_r)) * Y_reg_inv;
        else
            C_l(bin, :) = abs(hrtfs_l(bin, :)) .* exp(1i * phi_l) * Y_reg_inv;
            C_r(bin, :) = abs(hrtfs_r(bin, :)) .* exp(1i * phi_r) * Y_reg_inv;
        end

    % else f(bin) > fs/2*0.85
    %       do nothing

    end

end

% nyquist bin
C_l(end, :) = abs(C_l(end, :)); 
C_r(end, :) = abs(C_r(end, :));

% Do something with DC bin?

c_l_irs = ifft([C_l; conj(flipud(C_l(2:end-1, :)))], [], 1, 'symmetric');
c_r_irs = ifft([C_r; conj(flipud(C_r(2:end-1, :)))], [], 1, 'symmetric');

c_l_irs = fractional_delay(c_l_irs, size(c_l_irs, 1)/2 + grp_delay_hrirs_l);
c_r_irs = fractional_delay(c_r_irs, size(c_r_irs, 1)/2 + grp_delay_hrirs_r);

% ----------------------- window the irs, just in case --------------------
win      = hann(round_up_to_even(500*fs/48000)); 
fade_in  = win(1:end/2);
fade_out = win(end/2+1:end);

c_l_irs(1:length(fade_in),          :, :) = c_l_irs(1:length(fade_in),          :, :) .* repmat(fade_in,  1, size(c_l_irs, 2), size(c_l_irs, 3));
c_l_irs(end-length(fade_out)+1:end, :, :) = c_l_irs(end-length(fade_out)+1:end, :, :) .* repmat(fade_out, 1, size(c_l_irs, 2), size(c_l_irs, 3));

c_r_irs(1:length(fade_in),          :, :) = c_r_irs(1:length(fade_in),          :, :) .* repmat(fade_in,  1, size(c_r_irs, 2), size(c_r_irs, 3));
c_r_irs(end-length(fade_out)+1:end, :, :) = c_r_irs(end-length(fade_out)+1:end, :, :) .* repmat(fade_out, 1, size(c_r_irs, 2), size(c_r_irs, 3));

end
