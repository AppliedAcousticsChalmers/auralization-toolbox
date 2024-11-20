% This script converts HRIRs into an order-limited spherical harmonic
% representation by computing the spherical harmonics using MagLS.
%
% For best results at low frequencies, make sure that line 137 in 
% eMagLS/getMagLsFilters.m reads 'fade_win = getFadeWindow(len, .4);'.

clear;

addpath('../dependencies/'); 

% This is the repository from https://github.com/thomasdeppisch/eMagLS.git.
% It is included as a submodule.
addpath('eMagLS/'); 
addpath('eMagLS/dependencies/');
addpath('eMagLS/dependencies/Spherical-Harmonic-Transform/');

store_it = 1; % store the result 

N = 1; % desired SH order

% -------------------------------------------------------------------------

% a little tweak to the the low frequencies prettier
warning('Make sure that line 137 in getMagLsFilters.m reads ''fade_win = getFadeWindow(len, .4);''.');

% some metadata to be stored
comment = 'MagLS HRIR SH coefficients of Neumann KU100 Cologne data set';

taps = 512; % ir length

out_file_name = sprintf('hrirs_ku100_magls_sh_N%d.mat', N);

% ------------------------ load the raw HRIRs -----------------------------

% download the HRTFs if they don't exist
hrir_path = 'HRIR_L2702.sofa';
download_hrtfs(hrir_path);

SOFAstart;

hrirs_sofa = SOFAload(hrir_path);

% hrirs
hrirs_left  = squeeze(hrirs_sofa.Data.IR(:, 1, :)).';
hrirs_right = squeeze(hrirs_sofa.Data.IR(:, 2, :)).';

% incidence angles
hrir_grid_azi_rad = (hrirs_sofa.SourcePosition(:, 1)/180 * pi).';
hrir_grid_col_rad = (pi/2 - hrirs_sofa.SourcePosition(:, 2)/180 * pi).';

fs = double(hrirs_sofa.Data.SamplingRate);

% ------------------ convert to spherical harmonics -----------------------

[h_l_ring, h_r_ring] = getMagLsFilters(hrirs_left, hrirs_right, hrir_grid_azi_rad.', hrir_grid_col_rad.', N, fs, taps, false);    

% store the SH coefficients of the HRIRs
if store_it
    save(sprintf(out_file_name, N), 'h_l_ring', 'h_r_ring', 'N', 'fs', 'comment');
end

