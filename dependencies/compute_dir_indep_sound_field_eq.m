function [eq_ir] = compute_dir_indep_sound_field_eq(s_nm, hrirs_magls_sh, hrirs_orig_sofa, fs, N, ear)
% computes minimum phase EQ filter for the specified ear
% ear: 1 (left), 2 (right)

% do some zero padding to align lengths ( minimum length: 4096)
taps = max(size(s_nm, 1), 4096);
s_nm = [s_nm; zeros(taps-size(s_nm, 1), size(s_nm, 2))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% only EQ ipsilateral incidence directions
if ear == 1
    head_orientations_rad = flip(-linspace(0, pi/2, 21));
else
    head_orientations_rad =       linspace(0, pi/2, 21);
end

% ----- compute the deviation between auralization and ground truth -------

% because of convolution (must be even)
brir_length = size(s_nm, 1) + size(hrirs_magls_sh.h_l_ring, 1);

% un-eq'd auralization
brtfs    = zeros(brir_length/2+1, length(head_orientations_rad));
% ground truth
brtfs_gt = zeros(brir_length/2+1, length(head_orientations_rad));

% loop over head orientations 
for ii = 1 : length(head_orientations_rad)
    
    head_orientation_rad = head_orientations_rad(ii);

    % brirs will have more samples than s_nm
    brirs = render_ambi_signals_binaurally_t(s_nm, hrirs_magls_sh, [], head_orientation_rad, N);
    
    brtfs_tmp    = fft(brirs(:, ear), [], 1);
    brtfs(:, ii) = brtfs_tmp(1:end/2+1, :);
    
    % --- ground truth ---
    idx = SOFAfind(hrirs_orig_sofa, -head_orientation_rad/pi*180, 0);
    
    % extract irs
    hrirs_gt  = squeeze(hrirs_orig_sofa.Data.IR(idx(1), :, :)).';
    
    % zero padding 
    hrirs_gt = [hrirs_gt; zeros(brir_length-size(hrirs_gt, 1), size(hrirs_gt, 2))];
    
    brtfs_gt_tmp    = fft(hrirs_gt(:, ear), [], 1);
    brtfs_gt(:, ii) = brtfs_gt_tmp(1:end/2+1, :);
      
end

% --------------------------- plot it all ---------------------------------
% h = figure;
% set(gcf, 'Color', [1 1 1]);
% set(gcf, 'Position', [700 100 1600 500])

f = linspace(0, fs/2, size(brtfs, 1)).';

% subplot(1, 6, 1);
% 
% surf(head_orientations_rad/pi*180, f+5*eps, 20*log10(abs(brtfs_gt)), 'LineWidth', 2);
% view(0, 90);
% 
% shading interp;
% set(gca, 'YScale', 'log');
% 
% clim([-30 30]);
% xlim([head_orientations_rad(1) head_orientations_rad(end)]/pi*180);
% ylim([50 20000]);
% 
% xlabel('Azimuth');
% ylabel('Frequency');
% title('Ground truth');
% 
% subplot(1, 6, 2);
% 
% surf(head_orientations_rad/pi*180, f+5*eps, 20*log10(abs(brtfs)), 'LineWidth', 2);
% view(0, 90);
% 
% shading interp;
% set(gca, 'YScale', 'log');
% 
% clim([-30 30]);
% xlim([head_orientations_rad(1) head_orientations_rad(end)]/pi*180);
% ylim([50 20000]);
% 
% xlabel('Azimuth');
% ylabel('Frequency');
% title('Auralization');

% --- difference ---
diff = abs(brtfs) ./ abs(brtfs_gt);

% --- clip ---
%diff(diff > 2) = 2;

% --- smooth ---
smoothing_frac = 3;
diff_smoothed = AKfractOctSmooth(diff, 'amp', fs, smoothing_frac);
% --------------

% figure;
% semilogx(f+5*eps, 20*log10(abs(diff)), 'LineWidth', 2);
% hold on;
% semilogx(f+5*eps, 20*log10(abs(diff_smoothed)), 'LineWidth', 2);
% hold off;
% grid on
% xlim([50 20000]);

% subplot(1, 6, 3);
% 
% surf(head_orientations_rad/pi*180, f+5*eps, 20*log10(abs(diff)), 'LineWidth', 2);
% view(0, 90);
% 
% shading interp;
% set(gca, 'YScale', 'log');
% 
% clim([-30 30]);
% xlim([head_orientations_rad(1) head_orientations_rad(end)]/pi*180);
% ylim([50 20000]);
% 
% xlabel('Azimuth');
% ylabel('Frequency');
% title('Difference');
% 
% subplot(1, 6, 4);
% 
% surf(head_orientations_rad/pi*180, f+5*eps, 20*log10(abs(diff_smoothed)), 'LineWidth', 2);
% view(0, 90);
% 
% shading interp;
% set(gca, 'YScale', 'log');
% 
% clim([-30 30]);
% xlim([head_orientations_rad(1) head_orientations_rad(end)]/pi*180);
% ylim([50 20000]);
% 
% xlabel('Azimuth');
% ylabel('Frequency');
% title('Difference smoothed');

% --------------------- compute the eq filter -----------------------------

% average deviation between auralization and ground truth
spectral_signature = mean(diff_smoothed, 2);

eq = 1./spectral_signature;
% normalize
eq = eq./median(abs(eq(2:10)));

% no eq below 500 Hz
eq(f < 500) = 1;

% make minimum phase
eq_ir = ifft([eq; flipud(eq(2:end-1, :))], [], 1, 'symmetric');
[~, eq_ir] = rceps(eq_ir);

% shorten so that only high frequencies will be affected (130 taps at 48
% kHz sampling frequency)
tap_max = 2 * ceil(round(fs/48000 * 130)/2); % make sure it's even
eq_ir(tap_max:end) = 0;
 
eq_ir = eq_ir(1:tap_max);

eq = fft(eq_ir);
eq = eq(1:end/2+1, :);

% plot EQ

% % -------------------------
% f = linspace(0, fs/2, length(eq)).';
% figure;
% semilogx(f+5*eps, 20*log10(abs(eq)), ':', 'LineWidth', 2);
% grid on;
% xlim([30 20000]);
% ylim([-20 20]);
% xlabel('Frequency (Hz)');
% % -------------------------

% --------------------------- plot EQ'd -----------------------------------
% figure(h);
% 
% subplot(1, 6, 5);
% 
% surf(head_orientations_rad/pi*180, f+5*eps, 20*log10(abs(brtfs .* eq)), 'LineWidth', 2);
% view(0, 90);
% 
% shading interp;
% set(gca, 'YScale', 'log');
% 
% clim([-30 30]);
% xlim([head_orientations_rad(1) head_orientations_rad(end)]/pi*180);
% ylim([50 20000]);
% 
% xlabel('Azimuth');
% ylabel('Frequency');
% title('Ed''d');
% 
% subplot(1, 6, 6);
% 
% surf(head_orientations_rad/pi*180, f+5*eps, 20*log10(abs(brtfs .* eq)) - 20*log10(abs(brtfs_gt)), 'LineWidth', 2);
% view(0, 90);
% 
% shading interp;
% set(gca, 'YScale', 'log');
% 
% clim([-30 30]);
% xlim([head_orientations_rad(1) head_orientations_rad(end)]/pi*180);
% ylim([50 20000]);
% 
% xlabel('Azimuth');
% ylabel('Frequency');
% title('Ed''d difference');

end
