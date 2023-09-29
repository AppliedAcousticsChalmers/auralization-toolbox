function [ out ] = sphbesselj( nu, z )
%function [ out ] = sphbesselj( nu, z )
% SPHBESSELJ spherical bessel function of first type of order nu and 
% argument z

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

out = zeros( size( z ) );

% if more than one "nu"
if ( ( size(nu) > 1 ) )
    
    if ( size(nu) ~= size(z) )
        error('The sizes of NU and Z do not match.');
    end
    
    % evaluate for z~=0 to avoid division by "0"
    out( z~=0 ) = sqrt( pi ./ ( 2 .* z( z~=0 ) ) ) .* ...
                               besselj( nu( z~=0 ) + 0.5, z( z~=0 ) );
% if only one "nu"
else
    % evaluate for z~=0 to avoid division by "0"
    out( z~=0 ) = sqrt( pi ./ ( 2 .* z( z~=0 ) ) ) .* ...
                               besselj( nu + 0.5, z( z~=0 ) );
end
                                
% evaluate the rest
out( nu == 0 & z==0 ) = 1;
out( nu ~= 0 & z==0 ) = 0;

% if ( nu==0 )
%     out( find( z==0 ) ) = 1;
% elseif ( nu~=0 )
%     out( find( z==0 ) ) = 0;
% end
