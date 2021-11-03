function [] = draw_coil_profileview(xc,yc,W,cs)
% draw_coil_profileview  Draws a profile view of a coil
%
% [] = draw_coil_topview(xc,yc,W,cs)
%
% Inputs:
% xc       = center coordinate for x
% yc       = center coordinate for y
% W        = structure of winding coil parameters
%  W.lcr   = length of core (m)
%  W.cl    = clearance on each side of core length (m)
%  W.wcr   = width of core (m)
%  W.cw    = clearance on each side of core width (m)
%  W.wcl   = width of coil cross section (m)
%  W.dcl   = depth of coil cross section (m)
%  W.N     = number of turns
%  W.acnd  = area of conductor (m^2)
%  W.lcnd  = length of conductor (m)
%  W.vcnd  = volume of conductor (m^3)
%  W.pf    = packing factor
%  W.wwn   = width of window {defined as rectangle which is at least
%            W.cwn from the coil} (m)
%  W.lwn   = length of window (m)
%  W.cwn   = clearance to window (m)
%  W.lex   = mean excess length of bundle (beyond perimeter of core)
%  W.vol   = volume of coil (including dead space)
%  W.gvld  = 1 if valid geometry 0<=C.gvld<1 if not
% cs       = color string descriptor
%
% Written by:
% S.D. Sudhoff 
% School of Electrical and Computer Engineering
% 1285 Electrical Engineering Building
% West Lafayette, IN 47907-1285
% Phone: 765-494-3246
% Fax: 765-494-0676
% E-mail: sudhoff@ecn.purdue.edu

% center section
x0=-W.wwn/2-W.cwn;
y0= -W.dcl/2;
x1= x0;
y1= y0+W.dcl;
x2=-x0;
y2= y1;
x3=x2;
y3= y1-W.dcl;
xp=[x0 x1 x2 x3 x0];
yp=[y0 y1 y2 y3 y0];
fill(xc+xp,yc+yp,cs,xc-xp,yc+yp,cs,'EdgeColor',cs);

% cross sections
x0= -W.wwn/2-W.cwn-W.wcl;
y0= -W.dcl/2;
x1= x0;
y1= y0+W.dcl;
x2= x1+W.wcl;
y2= y1;
x3= x2;
y3= y1-W.dcl;
xp=[x0 x1 x2 x3 x0];
yp=[y0 y1 y2 y3 y0];
fill(xc+xp,yc+yp,cs,xc-xp,yc+yp,cs,'EdgeColor','k');

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

