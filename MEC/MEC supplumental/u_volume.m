function [vtot,vlegs,vbase] = u_volume(U)
% u_volume  Computes volume calculations on a U-core.
%
% [vtot,vlegs,vbase] = u_volume(U)
%
% Inputs:
% U       = stucture of uu core inductor parameters
%  U.du   = depth of U-core base (m)          
%  U.wu   = width of U-core leg (m)           
%  U.ws   = width of slot (m)                
%  U.ds   = depth of slot (m)
%  U.l    = length of U-core into page
%            Dimensions are such that height of U = du+ds
%                                     width of  U = 2*wu+ws
%                                     thickness of U = l
%
% Outputs:
%
% vtot     = volume of magnetic material parameters in m^3
% vlegs    = volume of both legs in m^3
% vbase    = voluem of base in m^3
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

vbase=(U.ws+U.wu)*U.ds*U.l;
vlegs=(2*U.ds+U.du)*U.wu*U.l;
vtot=(2*U.wu*U.ds+(U.ws+2*U.wu)*U.du)*U.l;

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
%  License along with MEC Toolbox.  If not, see
%  <http://www.gnu.org/licenses/>.
