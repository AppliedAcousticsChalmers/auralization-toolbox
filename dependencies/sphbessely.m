function [ out ] = sphbessely( nu, z )
%function [ out ] = sphbessely( nu, z )
%
%SPHBESSELY spherical bessel function of second kind of order nu and 
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

out         = zeros( size( z ) ); 
% avoid division by 0
out( z~=0 ) = sqrt( pi ./ ( 2.* z( z~=0 ) ) ) .* bessely( nu + 0.5, z( z~=0 ) );
out( z==0 ) = NaN;
