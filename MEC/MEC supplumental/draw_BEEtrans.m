function [] = draw_BEEtrans(BE,PBW,SBW,cww,fn)
% draw_BEEtrans  Draws a BEE-core transformer
%
% [] = draw_BEEtrans(BE,PBW,SBW)
%
% Inputs:
% BE       = structure of core parameters
% PBW      = structure of primary winding parameters
% SBW      = structure of secondary winding parameters
% cww      = clearance from winding to winding
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
draw_BE_topview(0,0,BE,'c');
draw_coil_topview(0,0,PBW,'r');
hold off;
title('BEE Transformer - Top View');
xlabel('x-position, m');
ylabel('z-position, m');
V1=max(abs(axis));
axis([-V1 V1 -V1 V1]);

% profile view
figure(fn+1)
hold on;
draw_coil_profileview(0, cww/2,-1,PBW,'r');
draw_coil_profileview(0,-cww/2, 1,SBW,'m');
draw_BE_profileview(0,0, 1,BE,'c');
draw_BE_profileview(0,0,-1,BE,'c');
hold off;
title('BEE Transformer - Profile View');
xlabel('x-position, m');
ylabel('y-position, m');
V2=max(abs(axis));
% 
V=max(V1,V2);
figure(fn)
axis([-V V -V V]);
figure(fn+1)
axis([-V V -V V]);
hold off;

end

function [] = draw_BE_topview(xc,yc,BE,cs)

% draw_BE_topview  Draws a top view of a BE-core
%
% [] = draw_BE_topview(xc,yc,BE,cs)
%
% Inputs:
% xc       = center coordinate for x (m)
% yc       = center coordinate for y (m)
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

% draw base
x0 = -BE.re;
y0 = -BE.rcp;
x1 = x0;
y1 = y0+2*BE.rcp;
x2 = x1+2*BE.re;
y2 = y1;
x3 = x2;
y3 = y2-2*BE.rcp;
xp=[x0 x1 x2 x3 x0];
yp=[y0 y1 y2 y3 y0];
fill(xc+xp,yc+yp,cs,xc-xp,yc+yp,cs,'EdgeColor','k');

% outline outer post
x0=-BE.re+BE.wop;
y0=0;
x1=x0;
y0=-BE.rcp;
xp=[x0 x1];
yp=[y0 y1];
plot(xc+xp,xc+yp,'k');
plot(xc-xp,yc+yp,'k');

% outline inner post
theta=linspace(0,2*pi,360);
xp=BE.rcp*cos(theta);
yp=BE.rcp*sin(theta);
plot(xc+xp,yc+yp,'k');

% outline hole in post
xp=BE.rh*cos(theta);
yp=BE.rh*sin(theta);
fill(xc+xp,yc+yp,cs,xc-xp,yc+yp,'w','EdgeColor','k');

end

function [] = draw_BE_profileview(xc,yc,ys,BE,cs)

% draw_BE_profileview  Draws a profile view of a BE-core
%
% [] = draw_BE_profileview(xc,yc,ys,BE,cs)
%
% Inputs:
% xc       = center coordinate for x (m)
% yc       = center coordinate for y (m)
% ys       = y stretch factor (set to -1 to flip orientation);
% BE       = structure of core parameters
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

% legs
x0 = -BE.re;
y0 = 0;
x1 = x0;
y1 = y0-BE.hp-BE.hb;
x2 = x1+2*BE.re;
y2 = y1;
x3 = x2;
y3 = y2+BE.hp+BE.hb;
x4 = BE.rop;
y4 = y3;
x5 = x4;
y5 = y4-BE.hp;
x6 = BE.rcp;
y6 = y5;
x7 = x6;
y7 = y6+BE.hp;
x8 = x7-2*BE.rcp;
y8 = y7;
x9 = x8;
y9 = y8-BE.hp;
x10=-BE.rop;
y10= y9;
x11 = x10;
y11 = y10+BE.hp;
xp=[x0 x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x0];
yp=[y0 y1 y2 y3 y4 y5 y6 y7 y8 y9 y10 y11 y0]*ys;
fill(xc+xp,yc+yp,cs,xc-xp,yc+yp,cs,'EdgeColor','k');

x0=BE.rh;
y0=-BE.hp-BE.hb;
x1=BE.rh;
y1=0;
xp=[x0 x1];
yp=[y0 y1]*ys;
plot(xc+xp,yc+yp,'k--');
plot(xc-xp,yc+yp,'k--');

x0=BE.rop;
y0=0;
x1=x0;
y1=y0-BE.hp;
x2=BE.rcp;
y2=y1;
xp=[x0 x1 x2];
yp=[y0 y1 y2]*ys;
plot(xc+xp,yc+yp,'k--');
plot(xc-xp,yc+yp,'k--');

end


function [] = draw_coil_topview(xc,yc,BW,cs)

% draw_coil_topview  Draws a top view of a bobbin-coil
%
% [] = draw_coil_topview(xc,yc,BW,cs)
%
% Inputs:
% xc       = center coordinate for x (m)
% yc       = center coordinate for y (m)
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

% draw base
theta=linspace(0,2*pi-eps,360);
xp1=BW.ro*cos(theta);
yp1=BW.ro*sin(theta);
theta=linspace(2*pi-eps,0,360);
xp2=BW.ri*cos(theta);
yp2=BW.ri*sin(theta);
xp=[xp1 xp2];
yp=[yp1 yp2];
fill(xc+xp,yc+yp,cs,'EdgeColor',cs);
plot(xc+xp1,yc+yp1,'k');
plot(xc+xp2,yc+yp2,'k');

end

function [] = draw_coil_profileview(xc,yc,ys,BW,cs)

% draw_coil_topview  Draws a top view of a bobbin-coil
%
% [] = draw_coil_topview(xc,yc,BW,cs)
%
% Inputs:
% xc       = center coordinate for x (m)
% yc       = center coordinate for y (m)
% ys       = y value stretch factor
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

% draw base
x0=-BW.ro;
y0= 0;
x1= BW.ro;
y1= 0;
x2= x1;
y2= y1-BW.h;
x3= x2-2*BW.ro;
y3= y2;
xp=[x0 x1 x2 x3 x0];
yp=[y0 y1 y2 y3 y0]*ys;
fill(xc+xp,yc+yp,cs,'EdgeColor','k');

x0=-BW.ri;
y0=0;
x1=x0;
y0=-BW.h;
xp=[x0 x1];
yp=[y0 y1]*ys;
plot(xp+xc,yp+yc,'k-',-xp+xc,yp+yc,'k-');

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
