function [ Ynm ] = sphharm(n, m, beta, alpha, type)
% spherical harmonics of degree n and order m.
%
% [Ynm] = sphharm(n, m, beta, alpha, type);
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

if (nargin < 5)
    type = 'complex';
end

if (n < 0)
    error('Degree(n) must not be negative.');
end

if (n < abs(m))
    warning('Absolute value of order m must be less than or equal to the degree n.'); 
    Ynm = zeros(size(alpha));
    return;
end

Lnm = asslegendre(n, abs(m), cos(beta));

factor_1 = (2*n + 1) / (4*pi);
factor_2 = factorial(n - abs(m)) ./ factorial(n + abs(m));

% complex spherical harmonics as used througout the book
if (strcmp(type, 'complex'))
    
    Ynm = (-1).^m .* sqrt(factor_1 .* factor_2) .* Lnm .* exp(1i .* m .* alpha);

% complex without canceling Condon-Shortley phase; same as scipy
elseif (strcmp(type, 'complex_wo_cs'))
    
    Lnm = asslegendre(n, m, cos(beta));
    
    factor_2 = factorial(n - m) ./ factorial(n + m);
    
    Ynm = sqrt(factor_1 .* factor_2) .* Lnm .* exp(1i .* m .* alpha);

elseif(strcmp(type, 'real')) % Ambisonics!
        
    if (m ~= 0)
        factor_1 = 2 * factor_1;
    end
    
    if (m < 0)
        Ynm = (-1).^m .* sqrt(factor_1 .* factor_2) .* Lnm .* sin(abs(m) .* alpha);
    else % m >= 0
        Ynm = (-1).^m .* sqrt(factor_1 .* factor_2) .* Lnm .* cos(m .* alpha);
    end
    
else
    error('Unknown type.');    
    
end

end
