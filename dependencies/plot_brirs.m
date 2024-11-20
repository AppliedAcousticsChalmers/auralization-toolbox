brtfs_0 = fft(brirs_0, [], 1);
brtfs_0 = brtfs_0(1:end/2+1, :);

brtfs_90 = fft(brirs_90, [], 1);
brtfs_90 = brtfs_90(1:end/2+1, :);


% make sure that the ground truth raw HRTFs are there
hrir_path = 'hrtfs/HRIR_L2702.sofa';
download_hrtfs(hrir_path);

SOFAstart;

hrirs_sofa = SOFAload(hrir_path);

idx = SOFAfind(hrirs_sofa, [0 90]-head_orientation_azimuth_deg, [0 0]); % (0, 0) and (90, 0)

% extract irs
hrirs_gt_0  = squeeze(hrirs_sofa.Data.IR(idx(1), :, :)).';
hrirs_gt_90 = squeeze(hrirs_sofa.Data.IR(idx(2), :, :)).';

% do some zero padding to make the graphs smoother
taps_plotting = 4096;

hrirs_gt_0  = [hrirs_gt_0;  zeros(taps_plotting-size(hrirs_gt_0, 1),  size(hrirs_gt_0,  2))];
hrirs_gt_90 = [hrirs_gt_90; zeros(taps_plotting-size(hrirs_gt_90, 1), size(hrirs_gt_90, 2))];

brtfs_gt_0 = fft(hrirs_gt_0, [], 1);
brtfs_gt_0 = brtfs_gt_0(1:end/2+1, :);

brtfs_gt_90 = fft(hrirs_gt_90, [], 1);
brtfs_gt_90 = brtfs_gt_90(1:end/2+1, :);

% -------------------------------- plot -----------------------------------
f1 = linspace(0, fs/2, size(brtfs_0,    1)).';
f2 = linspace(0, fs/2, size(brtfs_gt_0, 1)).';

handle_f = figure;
set(gcf, 'Color', [1 1 1]);
set(gcf, 'Position', [700 100 800 500]);

subplot(2, 2, 1);

semilogx(f2+5*eps, 20*log10(abs(brtfs_gt_0(:, 1))), 'k', 'LineWidth', 2);

hold on;

semilogx(f1+5*eps, 20*log10(abs(brtfs_0(:, 1))), 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2);
plot([f_transition f_transition], [-20 10], ':', 'Color', [1 1 1]*.5, 'Linewidth', 2);

hold off;
grid on;
box on;
xlim([30 20000]);
ylim([-20 10]);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
legend('Ground truth', 'Auralization', 'Location', 'southwest');
title(sprintf('%d°, left ear', 0-round(head_orientation_azimuth_deg)));

subplot(2, 2, 2);

semilogx(f2+5*eps, 20*log10(abs(brtfs_gt_90(:, 1))), 'k', 'LineWidth', 2);

hold on;

semilogx(f1+5*eps, 20*log10(abs(brtfs_90(:, 1))), 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2);
plot([f_transition f_transition], [-20 10], ':', 'Color', [1 1 1]*.5, 'Linewidth', 2);

hold off;
grid on;
box on;
xlim([30 20000]);
ylim([-20 10]);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
legend('Ground truth', 'Auralization', 'Location', 'southwest');
title(sprintf('%d°, left ear', 90-round(head_orientation_azimuth_deg)));

subplot(2, 2, 3);

semilogx(f2+5*eps, 20*log10(abs(brtfs_gt_0(:, 2))), 'k', 'LineWidth', 2);

hold on;

semilogx(f1+5*eps, 20*log10(abs(brtfs_0(:, 2))), 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2);
plot([f_transition f_transition], [-20 10], ':', 'Color', [1 1 1]*.5, 'Linewidth', 2);

hold off;
grid on;
box on;
xlim([30 20000]);
ylim([-20 10]);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
legend('Ground truth', 'Auralization', 'Location', 'southwest');
title(sprintf('%d°, right ear', 0-round(head_orientation_azimuth_deg)));

subplot(2, 2, 4);

semilogx(f2+5*eps, 20*log10(abs(brtfs_gt_90(:, 2))), 'k', 'LineWidth', 2);

hold on;

semilogx(f1+5*eps, 20*log10(abs(brtfs_90(:, 2))), 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2);
plot([f_transition f_transition], [-20 10], ':', 'Color', [1 1 1]*.5, 'Linewidth', 2);

hold off;
grid on;
box on;
xlim([30 20000]);
ylim([-20 10]);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
legend('Ground truth', 'Auralization', 'Location', 'southwest');
title(sprintf('%d°, right ear', 90-round(head_orientation_azimuth_deg)));
