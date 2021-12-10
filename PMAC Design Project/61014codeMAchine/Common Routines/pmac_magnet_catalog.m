function [M]=pmac_magnet_catalog(n)
% pmac_electrical_parameters assigns the magnet parameters based on
%                            a code for a magnet type. Used in Chapter 9 
%                            of Power Magnetic Devices: A Multi-Objective 
%                            Design Approach," by S.D. Sudhoff
% Call:
% M=pmac_magnet_catalog(n)
%
% Inputs:
% n        = permanent magnet number
%            n=1  - NdFeB N35
%            n=2  - NdFeB N50
%            n=3  - NdFeB Plastic
%            n=4  - SmCo R20
%            n=5  - SmCo R32
%            n=6  - Ferrite AC-12
%            n=7  - AlNiCo 8H
%
% Outputs:
% M        = structure of permanent magnet parameters
%  M.desc  = magnet desription
%  M.br    = residual flux density (T)
%  M.x     = susceptibility 
%  M.Hlim  = minimum field intensity before demag (A/m)
%            (Taken to be 1/2 the intrinsic coercive force at room temp)
%  M.row   = mass density (kg/m^3)
%
% Written by:
% S.D. Sudhoff
% School of Electrical and Computer Engineering
% 1285 Electrical Engineering Building
% West Lafayette, IN 47907-1285
% Phone: 765-494-3246
% Fax: 765-494-0676
% E-mail: sudhoff@ecn.purdue.edu

switch n
   case 1
      M.desc='NdFeB N35';
      M.br  =1.19;
      M.Hci =-867e3;
      M.x   =0.09;
      M.row =7500;
   case 2
      M.desc='NdFeB N50';
      M.br  = 1.43;
      M.Hci = -836e3;
      M.x   = 0.36;
      M.row = 7500;
   case 3
      M.desc='NdFeB Plastic';
      M.br  = 0.66;
      M.Hci = -577e3;
      M.x   = 0.24;
      M.row = 5700;
   case 4
      M.desc='SmCo R20';
      M.br  = 0.9;
      M.Hci = -2400e3;
      M.x   = 0.02;
      M.row = 8400;
   case 5
      M.desc='SmCo R32';
      M.br  = 1.15;
      M.Hci = -1350e3;
      M.x   = 0.10;
      M.row = 8300;
   case 6
      M.desc='Ferrite AC-12';
      M.br  = 0.4;
      M.Hci = -318e3;
      M.x   = 0.10;
      M.row = 4900;
   case 7
      M.desc='AlNiCo 8H';
      M.br  = 0.74;
      M.Hci = -151e3;
      M.x   = 1.5;
      M.row = 7250;
   otherwise
      error('Unknown magnet type');
end
          
end