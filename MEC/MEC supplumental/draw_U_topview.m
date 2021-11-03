function [] = draw_U_topview(xc,yc,U,cs)

% draw_U_topview  Draws a top view of a U-core
%
% [] = draw_U_topview(xc,yc,U,cs)
%
% Inputs:
% xc       = center coordinate for x (m)
% yc       = center coordinate for y (m)
% U        = structure of core parameters
%  U.du    = depth of U-core base (m)          
%  U.wu    = width of U-core leg (m)           
%  U.ws    = width of slot (m)                
%  U.ds    = depth of slot (m)
%  U.l     = length of U-core into page
%            Dimensions are such that height of U = du+ds
%                                     width of  U = 2*wu+ws
%                                     thickness of U = l
%  U.pmub  = pointer to function which returns absolute permeabilty as 
%            function of flux density
%  U.MMP   = structure of core magnetic material parameters

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

% center
x0 = -U.ws/2;
y0 = -U.l/2;
x1 = x0;
y1 = y0+U.l;
x2 = x1+U.ws;
y2 = y1;
x3 = x2;
y3 = y2-U.l;
xp=[x0 x1 x2 x3 x0];
yp=[y0 y1 y2 y3 y0];
fill(xc+xp,yc+yp,cs,xc-xp,yc+yp,cs,'EdgeColor',cs);

% left side
x0 = -U.ws/2;
y0 = -U.l/2;
x1 = x0-U.wu;
y1 = y0;
x2 = x1;
y2 = y1+U.l;
x3 = x2+U.wu;
y3 = y2;
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
