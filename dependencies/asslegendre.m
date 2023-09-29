function [Lnm] = asslegendre(n, m, arg)
%ASSLEGENDRE Calculates associated Legendre function of degree n and order
% m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This work is supplementary material for the book                        %
%                                                                         %
% Jens Ahrens, Analytic Methods of Sound Field Synthesis, Springer-Verlag %
% Berlin Heidelberg, 2012, http://dx.doi.org/10.1007/978-3-642-25743-8    %
%                                                                         %
% It has been downloaded from http://soundfieldsynthesis.org and is       %
% licensed under a Creative Commons Attribution-NonCommercial-ShareAlike  %
% 3.0 Unported License. Please cite the book appropriately if you use     %
% these materials in your own work.                                       %
%                                                                         %
% (c) 2012 by Jens Ahrens                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (n < 0)
    
    % Gumerov and Duraiswami, Eq. (2.1.44)
    Lnm = asslegendre(abs(n)-1, m, arg);
    
    return;
end

if (n < abs(m))
    % warning('Absolute value of order(m) must be less than or equal to the degree(n).'); 
    Lnm = zeros(size(arg));
    return;
end

Lnm = legendre(n, arg);

if (n ~= 0)
    Lnm = squeeze(Lnm(abs(m) + 1, :, :));  
end

if (m < 0)
    Lnm = (-1).^abs(m) .* factorial(n - abs(m)) ./ factorial(n + abs(m)) .* Lnm;
end

Lnm = reshape(Lnm, size(arg));

end