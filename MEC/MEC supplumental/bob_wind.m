function [BW] = bob_wind(ri,ro,h,N,ac,aw,CP)
% bob_wind  Sets up data structure associated with a bobbin winding.
%
% [BW] = bob_wind(ri,ro,h,N,ac,aw)
% [BW] = bob_wind(ri,ro,h,N,ac,aw,CP);
%
% Inputs:
%  ri     = inner radius of bobbin (m)
%  ro     = outer radius of bobbin (m)
%  h      = height of bobbin (m)
%  N      = number of turns
%  ac     = conductor area (m^2)
%  aw     = wire area (m^2)
%  CP     = structure of conductor parameters
%
% Outputs:
% BW        = structure of bobin winding parameters
%  BW.ri    = inner radius of bobbin (m)
%  BW.ro    = outer radius of bobbin (m)
%  BW.w     = width of bobbin cross section (m)
%  BW.l     = mean length of bobbin (circumference of center) (m)
%  BW.h     = height of bobbin (m)
%  BW.N     = number of turns
%  BW.ac    = conductor area (m^2)
%  BW.aw    = wire area (m^2)
%  BW.bcs   = bobbin cross section (m^2)
%  BW.wpf   = wire packing factor 
%  BW.cpf   = conductor packing factor
%  BW.vc    = volume of conductor (m^3)
%  BW.lc    = length of conductor (m^3)  
%  BW.mc    = mass of conductor (m^3) {if CP supplied }
%  BW.rdc   = dc resistance of conductor (Ohms) {if CP supplied}
%
% Internal:
%
% Written by:
% S.D. Sudhoff
% School of Electrical and Computer Engineering
% 1285 Electrical Engineering Building
% West Lafayette, IN 47907-1285
% Phone: 765-494-3246
% Fax: 765-494-0676
% E-mail: sudhoff@ecn.purdue.edu

% assign known parameters
BW.ri=ri;
BW.ro=ro;
BW.h=h;
BW.N=N;
BW.ac=ac;
BW.aw=aw;

% compute bobbin cross section and cross section width
BW.l=pi*(BW.ro+BW.ri);
BW.w=BW.ro-BW.ri;
BW.bcs=BW.w*BW.h;

% compute packing factors
BW.wpf=BW.aw*BW.N/BW.bcs;
BW.cpf=BW.ac*BW.N/BW.bcs;

% compute conductor volume and length
BW.vc=pi*(BW.ro^2-BW.ri^2)*BW.h*BW.cpf;
BW.lc=BW.vc/BW.ac;

% optional calculations
if nargin==7
   BW.CP=CP;
   BW.mc=BW.vc*CP.row;
   BW.rdc=BW.lc/(BW.ac*CP.sigmac);
end

%  Copyright 2013 - Scott Sudhoff 
% 
%  The program is distributed under the terms of the GNU Lesser
%  General Public License(GNU LGPL). 
% 
%  This file is part of MEC Toolbox 3.2
% 
%  MEC Toolbox is free software: you can redistribute it and/or modify
%  it under the terms of the GNU Lesser General Public License as published
%  by the Free Software Foundation, either version 3 of the License, or 
%  (at your option) any later version.
% 
%  MEC Toolbox is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
%  See the GNU Lesser General Public License for more details.
% 
%  You should have received a copy of the GNU Lesser General Public
%  License along with MEC Toolbox.  If not, see <http://www.gnu.org/licenses/>.
