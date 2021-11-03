function [] = draw_UUtrans(U,PW,SW,fn)
% draw_UUtrans  Draws a UU-core transformer
%
% [] = draw_UUtrans(U,PW,SW)
%
% Inputs:
% U        = structure of core parameters
% PW       = structure of primary winding parameters
% SW       = structure of secondary winding parameters
% fn       = figure number
%
% Internal:
% x_pcc    = x-coordinate of center of primary winding
% x_scc    = x-coordinate of center of secondary winding
% V1,V2,V  = used to set axis limits
%
% Written by:
% S.D. Sudhoff
% School of Electrical and Computer Engineering
% 1285 Electrical Engineering Building
% West Lafayette, IN 47907-1285
% Phone: 765-494-3246
% Fax: 765-494-0676
% E-mail: sudhoff@ecn.purdue.edu


% top view
figure(fn);
hold on;
draw_U_topview(0,0,U,'c');
x_pcc=-U.ws/2-U.wu/2;
x_scc= U.ws/2+U.wu/2;
draw_coil_topview(x_pcc,0,PW,'r');
draw_coil_topview(x_scc,0,SW,'m');
hold off;
title('UU Transformer - Top View');
xlabel('x-position, m');
ylabel('z-position, m');
V1=max(abs(axis));

% profile view
figure(fn+1)
hold on;
draw_coil_profileview(x_pcc,0,PW,'r')
draw_coil_profileview(x_scc,0,SW,'m')
draw_U_profileview(0,0, 1,U,'c');
draw_U_profileview(0,0,-1,U,'c');
hold off;
title('UU Transformer - Profile View');
xlabel('x-position, m');
ylabel('y-position, m');
V2=max(abs(axis));

V=max(V1,V2);
figure(fn)
axis([-V V -V V]);
figure(fn+1)
axis([-V V -V V]);

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

