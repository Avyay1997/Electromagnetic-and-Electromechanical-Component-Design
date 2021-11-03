function [Pint] = Pint_rct_bndle_air(lcs1,lcs2,l)
% Pint_rct_bndle_air   Calculates the internal permeance of a rectangular
%                      wire bundle in air using method set forth in
%                      Brandon Cassimere, "Modeling and Evolutionary Design
%                      of Permanent Magnet Synchronous Machines," Phd
%                      Dissertation, Purdue University, May 2008.
%
% Pint = Pint_rct_bndle_air(lcs1,lcs2,l,N)
%
% Inputs:
% lcs1 = length of one side of the rectuangular cross section (m)
% lcs2 = length of other side of the rectangular cross section (m)
% l    = axial length of wire bundle
%
% Outputs:
% Pint = internal permeance associated with the bundle (H^-1)
%
% Written by:
% S.D. Sudhoff
% School of Electrical and Computer Engineering
% 1285 Electrical Engineering Building
% West Lafayette, IN 47907-1285
% Phone: 765-494-3246
% Fax: 765-494-0676
% E-mail: sudhoff@ecn.purdue.edu

mu0=4e-7*pi;

if (lcs1>lcs2)
   lr=lcs1;
   dr=lcs2;
else
   lr=lcs2;
   dr=lcs1;
end

if (lr==dr)
   Pint=mu0*l/32;
else
   dr2=dr*dr;
   dr3=dr2*dr;
   dr4=dr3*dr;
   lrmdr=lr-dr;
   lrmdr2=lrmdr*lrmdr;
   lrmdr3=lrmdr2*lrmdr;
   lrmdr4=lrmdr3*lrmdr;
   dr2lr2=dr2*lr^2;
   Pint=mu0*l*(dr4/32.0 + ...
               dr3*lrmdr/16.0 + ...
               dr2*lrmdr2/64.0 - ...
               dr*lrmdr3/64.0 + ...
               lrmdr4*log(1+2*dr/lrmdr)/128.0)/dr2lr2;
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
%  License along with MEC Toolbox.  If not, see
%  <http://www.gnu.org/licenses/>.
