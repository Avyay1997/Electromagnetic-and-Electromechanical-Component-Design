function [C] = rect_coil(lcr,cl,wcr,cw,dcl,wcl,Nt,Np,ac)
% rect_coil  Calculates dimensions of a rectangular core
%
% [C] = rect_coil(lcr,cl,wcr,cw,dcl,wcl,Nt,Np,ac)
%
% Inputs:
%  lcr     = length of core segment the coil goes around (m)
%  cl      = clearance on each side of the length (m) 
%  wcr     = width of core segment the coil goes around (m)
%  cw      = clearance on each side of the width (m)
%  dcl     = depth of coil cross section (m)
%  wcl     = width of coil cross section (m)
%  Nt      = number of turns 
%  Np      = number of parallel conductor paths in making turns
%  ac      = area of conductor (m^2)
%
% Outputs:
% C        = structure of coil parameters
%  C.lcr   = length of core (m)
%  C.cl    = clearance on each side of core length (m)
%  C.wcr   = width of core (m)
%  C.cw    = clearance on each side of core width (m)
%  C.wcl   = width of coil cross section (m)
%  C.dcl   = depth of coil cross section (m)
%  C.Nt    = number of turns
%  C.Np    = number of parallel paths
%  C.acnd  = area of conductor (m^2)
%  C.lcnd  = total length of conductor (m)
%            {all the paralllel paths stacked on end}
%  C.vcnd  = volume of conductor (m^3)
%  C.pf    = packing factor
%  C.wwn   = width of window {defined as rectangle which is at least
%            C.cwn from the coil} (m)
%  C.lwn   = length of window (m)
%  C.cwn   = clearance to window (m)
%  C.lex   = mean excess length of bundle (beyond perimeter of core)
%  C.vol   = volume of coil (including dead space)
%  C.gvld  = 1 if valid geometry 0<=C.gvld<1 if not
%
% Written by:
% S.D. Sudhoff
% School of Electrical and Computer Engineering
% 1285 Electrical Engineering Building
% West Lafayette, IN 47907-1285
% Phone: 765-494-3246
% Fax: 765-494-0676
% E-mail: sudhoff@ecn.purdue.edu

% assign known information to fields
C.lcr=lcr;
C.cl=cl;
C.wcr=wcr;
C.cw=cw;
C.Nt=Nt;
C.Np=Np;
C.acnd=ac;
C.wcl=wcl;
C.dcl=dcl;

% compute the bending radius and winding window
C.cwn=max(cw,cl);
C.wwn=wcr+2*cw-2*C.cwn;
C.lwn=lcr+2*cl-2*C.cwn;
if (C.wwn>0)&&(C.lwn>0)
   C.gvld=1;
else
   C.gvld=exp(min(C.wwn,C.lwn))-1.0;
end

% compute the mean excess lenth of coil
C.lex=2*C.lwn+2*C.wwn+2*pi*(C.cwn+0.5*C.wcl)-2*lcr-2*wcr;

% compute the volume of the coil
C.vol=(2.0*C.wcl*C.wwn+...
       2.0*C.wcl*C.lwn+...
       pi*((C.wcl+C.cwn).^2-C.cwn.^2))*C.dcl;
   
% compute the length of the conductor
C.pf=Nt*Np*ac/(C.wcl*C.dcl);
C.vcnd=C.vol*C.pf;
C.lcnd=C.vcnd/ac;

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
