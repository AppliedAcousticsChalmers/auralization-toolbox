function [ out ] = sphbesselh( nu, k, z ) 
% function [ out ] = sphbesselh( nu, k, z )
%
% SPHBESSELH spherical hankel function of order nu, kind k, and argument z

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

if ( k==1 ) 
    sign = 1;
elseif ( k==2 )
    sign = -1;
else
    error('Invalid kind of Hankel function is asked.');
end

out = sphbesselj( nu, z ) + 1i .* sign .* sphbessely( nu, z );
