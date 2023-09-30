function [C, condition_number] = get_c_nm_direct(irs_calibration, hrirs_calibration, fs, fc_reg, reg_parameter)

fc_bp = [10 18000];

tfs = fft(irs_calibration, [], 1);
tfs = tfs(1:end/2+1, :, :);

% zero pad HRIRs
hrirs_calibration = [hrirs_calibration; ...
    zeros(size(irs_calibration, 1) - size(hrirs_calibration, 1), size(hrirs_calibration, 2))];

hrtfs = fft(hrirs_calibration, [], 1);
hrtfs = hrtfs(1:end/2+1, :);

condition_number = zeros(size(tfs, 1), 1);
C = zeros(size(tfs, 1), size(tfs, 2));

f = linspace(0, fs/2, size(tfs, 1)).';

for bin = 1 : size(tfs, 1)
    
    if f(bin) < fc_reg(1)
        lambda = reg_parameter(2);
    elseif f(bin) > fc_reg(2)
        lambda = reg_parameter(2);
    else 
        lambda = reg_parameter(1);
    end
    
    % source directions x mics
    V = squeeze(tfs(bin, : , :)).';
    
    % hrtf
    h = hrtfs(bin, :).';
    
    matrix = V'*V + lambda * eye(size(V, 2));
    
    C(bin, :) = (inv(matrix) * V' * h).';
    
    condition_number(bin) = cond(matrix);
        
end

c_irs = ifft([C; conj(flipud(C(2:end-1, :)))], [], 1, 'symmetric');

% ------------------------- bandpass filtering ----------------------------


% highpass
[z,p,k] = butter(9, fc_bp(1)/(fs/2), 'high');
sos_hp = zp2sos(z,p,k);

%fvtool(sos_hp,'Analysis','freq'); xlim([0 0.02]);
%zplane(z, p);
%grid on;

% lowpass
[z,p,k] = butter(9, fc_bp(2)/(fs/2));
sos_lp = zp2sos(z,p,k);
%fvtool(sos_lp,'Analysis','freq');

c_irs = sosfilt(sos_hp, c_irs);
c_irs = sosfilt(sos_lp, c_irs);

C = fft(c_irs, [], 1);
C = C(1:end/2+1, :);

end