function [c_l_irs, c_r_irs, condition_number] = get_c_direct(irs_calibration, hrirs_calibration_l, hrirs_calibration_r, fs, f_transition, dynamic_range_permitted_dB)
%
% irs_calibration: input signal
% hrirs_calibration_l, hrirs_calibration_r: desired output signal

fprintf('Preparing computation of the auralization matrices ... ');

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

fprintf('done.\n\n');

condition_number = zeros(size(tfs, 1), 1);

C_l = zeros(size(tfs, 1), size(tfs, 2));
C_r = zeros(size(tfs, 1), size(tfs, 2));

% -------------------------------------------------------------------------
display_progress(sprintf('Computing the auralization matrices (using eMagLS2 above %d Hz):', round(f_transition)));

% loop over all bins
for bin = 1 : size(tfs, 1)

    display_progress(bin/size(tfs, 1));

    pw_grid = squeeze(tfs(bin, :, :));

    %[U, s, V] = svd(pw_grid.', 'econ', 'vector');
    [U, s, V] = svd(pw_grid.', 'econ'); s = diag(s); % for older MATLAB versions

    dynamic_range_actual_dB   = 20*log10(s(1)/s(end));
    dynamic_range_to_apply_dB = get_dynamic_range_to_apply_dB(f(bin), dynamic_range_actual_dB, dynamic_range_permitted_dB);

    % reduce rank by 1 for safety
    s = s(1:end-1);
    U = U(:, 1:end-1);
    V = V(:, 1:end-1);

    % invert the matrix
    s = 1 ./ max(s, 10^(-dynamic_range_to_apply_dB/20) * max(s));
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

fprintf('\n\n');

end

% -------------------------------------------------------------------------
function dynamic_range_to_apply_dB = get_dynamic_range_to_apply_dB(f, dynamic_range_actual_dB, dynamic_range_permitted_dB)

head_room_dB = 3; % make sure that some amount of regularization is always being applied

% --- determine how much to regularize ---
if f > 10000 % 10 kHz

    if dynamic_range_actual_dB < dynamic_range_permitted_dB(3) + head_room_dB
        dynamic_range_to_apply_dB = dynamic_range_actual_dB - head_room_dB;
    else
        dynamic_range_to_apply_dB = dynamic_range_permitted_dB(3);
    end

elseif f > 100 % 100 Hz - 10 kHz

    if dynamic_range_actual_dB < dynamic_range_permitted_dB(2) + head_room_dB
        dynamic_range_to_apply_dB = dynamic_range_actual_dB - head_room_dB;
    else
        dynamic_range_to_apply_dB = dynamic_range_permitted_dB(2);
    end

else
    
    if dynamic_range_actual_dB < dynamic_range_permitted_dB(1) + head_room_dB
        dynamic_range_to_apply_dB = dynamic_range_actual_dB - head_room_dB;
    else
        dynamic_range_to_apply_dB = dynamic_range_permitted_dB(1);
    end

end

dynamic_range_to_apply_dB = max(0, dynamic_range_to_apply_dB);

end
% -------------------------------------------------------------------------