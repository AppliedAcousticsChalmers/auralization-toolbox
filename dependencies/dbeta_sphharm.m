function [Ynm] = dbeta_sphharm(n, m, beta, alpha)
% Derivative of spherical harmonic of degree n and order m with respect to
% the colatitude. Only available for the type 'real'.
%
% [Ynm] = dalpha_sphharm(n, m, beta, alpha);
%
% n           - spherical harmonic degree
% m           - spherical harmonic order
% beta        - colatitude to be calculated
% alpha       - azimuth to be calculated 
% type        - 'complex' (default), 'complex_wo_cs', ' 'real'
%
% alpha and beta can be arrays but have to be of same size or one of them
% has to be a scalar.
%
% Written by Jens Ahrens, 2022

if (n < 0)
    error('Degree(n) must not be negative.');
end

if (n < abs(m))
    warning('Absolute value of degree m must be less than or equal to the order n.'); 
    Ynm = zeros(size(alpha));
    return;
end

Lnm_prime = 1./(1-cos(beta).^2) .* ((n+1)*(n+abs(m))/(2*n+1) .* asslegendre(n-1, abs(m), cos(beta)) ...
                                    - n*(n-abs(m)+1)/(2*n+1) .* asslegendre(n+1, abs(m), cos(beta)));

Lnm_prime((1-cos(beta).^2) == 0) = 0;                               
                                
                                
factor_1 = (2*n + 1) / (4*pi);
factor_2 = factorial(n - abs(m)) ./ factorial(n + abs(m));
    
if (m ~= 0)
    factor_1 = factor_1 * 2; % This is the factor sqrt(2)
end

if (m < 0)
    Ynm = (-1).^m .* sqrt(factor_1 .* factor_2) .* (-1) .* sin(beta) .* Lnm_prime .* sin(abs(m) .* alpha);
else % m >= 0
    Ynm = (-1).^m .* sqrt(factor_1 .* factor_2) .* (-1) .* sin(beta) .* Lnm_prime .* cos(abs(m) .* alpha);
end


end
