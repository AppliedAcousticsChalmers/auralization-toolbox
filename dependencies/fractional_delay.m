function [sig] = fractional_delay(sig, delay_samples)
% applies a fractional delay in the frequency domain

% normalized angular frequency

omega   = linspace(0, 0.5, size(sig, 1)/2+1).';
delay_f = exp(-1i * 2*pi * omega .* delay_samples);

% enforce Nyquist bin to be real
delay_f(end, :, :) = real(delay_f(end, :, :)); 

% double sided spectrum
delay_f = [delay_f; flipud(conj(delay_f(2:end-1, :, :)))];

% apply complex delay in frequency domain
sig_f = fft(sig);
sig_f = sig_f .* delay_f;
sig   = ifft(sig_f, [], 1, 'symmetric');

end