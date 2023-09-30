function [sh_coeffs_t] = get_sh_coefficients_t(irs, c_nm, N)

fprintf('Computing the SH coefficients ... ');

sh_coeffs_t = zeros(size(irs, 1), (N+1)^2);

irs = double(irs); % fftfilt requires this

% loop over all SH modes
for n = 0 : N
   for m = -n : n 
       
       % apply the quadrature matrix to the simulation data
       sh_coeffs_t(:, n^2+n+m+1) = sum(fftfilt(double(c_nm(:, :, n^2+n+m+1)), irs), 2);
       
   end
end

fprintf('done.\n\n');

end

