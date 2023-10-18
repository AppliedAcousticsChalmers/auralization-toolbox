function [out] = d_sh_mode(n, m, r, beta, alpha, k, type)
% Derivative of j_n(kr) * Y_n,m(beta, alpha) w.r.t. the three Cartesian
% coordinate dimensions x, y, and z
%
% type: 'dx', 'dy', 'dz'

if (strcmp(type, 'dx'))
    
    factor_r     =   cos(alpha) .* sin(beta);
    factor_beta  =   cos(alpha) .* cos(beta);
    factor_alpha = - sin(alpha);
    
elseif (strcmp(type, 'dy'))
    
    factor_r     =   sin(alpha) .* sin(beta);
    factor_beta  =   sin(alpha) .* cos(beta);
    factor_alpha =   cos(alpha);
 
elseif (strcmp(type, 'dz'))
    
    factor_r     =   cos(beta);
    factor_beta  = - sin(beta);
    factor_alpha =   zeros(size(beta));
end

bessel       = sphbesselj(n, k.*r);
bessel_prime = 1/(2*n+1) * (n * sphbesselj(n-1, k.*r) - (n+1) * sphbesselj(n+1, k.*r));

% dr (this includes the poles)
out = factor_r .* ...
      k .* bessel_prime .* sphharm(n, m, beta, alpha, 'real');

% --- compute the following only if we do not evalute the derivative at a pole ---
indices = beta ~= 0 & beta ~= pi;

% dbeta
out(indices) = out(indices) + factor_beta(indices) .* ...
      1./r(indices) .* bessel(indices) .* dbeta_sphharm(n, m, beta(indices), alpha(indices));
  
% dalpha
out(indices) = out(indices) + factor_alpha(indices) .* ...
      1./(r(indices) .* sin(beta(indices))) .* bessel(indices) .* dalpha_sphharm(n, m, beta(indices), alpha(indices));

% edge case
if find(r == 0)
    warning('I cannot evaluate the spatial gradient in the coordinate origin.');
    out(r == 0) = 0;
end
  
end