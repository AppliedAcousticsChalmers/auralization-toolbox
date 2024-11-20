function [sh_coeffs_t] = get_sh_coefficients_t(irs, c_nm, N)

irs = double(irs); % fftfilt requires this

% zero pad for the convolution
irs = [irs; zeros(size(c_nm, 1), size(irs, 2))];

sh_coeffs_t = zeros(size(irs, 1), (N+1)^2);

% loop over all SH modes
for n = 0 : N
   for m = -n : n 
       
       % apply the quadrature matrix to the simulation data
       sh_coeffs_t(:, n^2+n+m+1) = sum(fftfilt(double(c_nm(:, :, n^2+n+m+1)), irs), 2);
       
   end
end

end

