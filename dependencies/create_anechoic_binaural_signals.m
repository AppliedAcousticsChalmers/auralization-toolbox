% get an example audio signal
[sig, fs_sig] = audioread('resources/drum_loop_48k.wav');

assert(fs == fs_sig);

% render audio signal binaurally
out_brirs_0     = fftfilt(brirs_0,     sig);
out_brirs_90    = fftfilt(brirs_90,    sig);
out_brirs_gt_0  = fftfilt(hrirs_gt_0,  sig);
out_brirs_gt_90 = fftfilt(hrirs_gt_90, sig);

max_amp = max(abs([out_brirs_0(:); out_brirs_90(:); out_brirs_gt_0(:); out_brirs_gt_90(:)]));

% normalize
out_brirs_0     = out_brirs_0     ./ max_amp;
out_brirs_90    = out_brirs_90    ./ max_amp;
out_brirs_gt_0  = out_brirs_gt_0  ./ max_amp;
out_brirs_gt_90 = out_brirs_gt_90 ./ max_amp;

% store

fprintf('Storing anechoic audio examples in the files\n')
fprintf('   auralization_binaural_preview_anechoic_0deg.wav\n')
fprintf('   auralization_binaural_preview_anechoic_90deg.wav\n')
fprintf('   ground_truth_binaural_anechoic_0deg.wav\n')
fprintf('   ground_truth_binaural_anechoic_90deg.wav\n\n')

audiowrite('auralization_binaural_preview_anechoic_0deg.wav',  out_brirs_0, fs);
audiowrite('auralization_binaural_preview_anechoic_90deg.wav', out_brirs_90, fs);
audiowrite('ground_truth_binaural_anechoic_0deg.wav',  out_brirs_gt_0, fs);
audiowrite('ground_truth_binaural_anechoic_90deg.wav', out_brirs_gt_90, fs);

fprintf('Done.\n\n');

