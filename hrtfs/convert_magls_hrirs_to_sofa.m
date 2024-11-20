clear;

addpath('../dependencies/');

% download SOFA template file if it does not exist
hrir_path = 'HRIR_L2702.sofa';
download_hrtfs(hrir_path);

SOFAstart;

hrirs_sofa = SOFAload(hrir_path);

% store to later comparison
hrirs_sofa_ref = hrirs_sofa;

% incidence angles
azi_rad = (hrirs_sofa.SourcePosition(:, 1)/180 * pi).';
col_rad = (pi/2 - hrirs_sofa.SourcePosition(:, 2)/180 * pi).';

% loop over all HRTF sets
for N = 1 : 25

    hrirs_sh = load(sprintf('hrirs_ku100_magls_sh_N%d.mat', N));

    hrirs_left  = zeros(size(hrirs_sh.h_l_ring, 1), length(hrirs_sofa.SourcePosition(:, 1)));
    hrirs_right = zeros(size(hrirs_left));

    for n = 0 : N
        for m = -n : n

            hrirs_left  = hrirs_left  + hrirs_sh.h_l_ring(:, n^2+n+m+1) .* sphharm(n, m, col_rad, azi_rad, 'real');
            hrirs_right = hrirs_right + hrirs_sh.h_r_ring(:, n^2+n+m+1) .* sphharm(n, m, col_rad, azi_rad, 'real');                

        end
    end

    % convert to SOFA
    hrirs_sofa.Data.IR = single(cat(2, permute(hrirs_left, [2 3 1]), permute(hrirs_right, [2 3 1])));

    hrirs_sofa = SOFAupdateDimensions(hrirs_sofa, 'verbose', 1);

    SOFAsave(sprintf('hrirs_ku100_magls_N%d.sofa', N), hrirs_sofa);

end

% test a few IRs
taps  = 4096;
index = 1000;

hrirs_ref  = squeeze(hrirs_sofa_ref.Data.IR(index, 1, :));
hrirs_test = squeeze(hrirs_sofa.Data.IR(    index, 1, :));

% zero pad
hrirs_ref  = [hrirs_ref;  zeros(taps-size(hrirs_ref,  1), size(hrirs_ref,  2))];
hrirs_test = [hrirs_test; zeros(taps-size(hrirs_test, 1), size(hrirs_test, 2))];

fs = double(hrirs_sofa_ref.Data.SamplingRate);

f_ref  = linspace(5*eps, fs, size(hrirs_ref , 1)).';
f_test = linspace(5*eps, fs, size(hrirs_test, 1)).';

figure;
semilogx(f_ref,  20*log10(abs(fft(hrirs_ref ))), 'LineWidth', 2);
hold on;
semilogx(f_test, 20*log10(abs(fft(hrirs_test))), ':', 'LineWidth', 2);
hold off;

grid on;
xlim([30 20000]);
ylim([-30 20]);


