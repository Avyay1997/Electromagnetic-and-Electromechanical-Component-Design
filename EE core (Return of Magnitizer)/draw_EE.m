function [] = draw_EE(EE,W,ksc,fn)
% draw_EE draws EI core inductor
%
% Call:
% draw_EE(EE,W,ksc,fn)
%
% Inputs:
% EE       = structure of EE core parameters
%  EE.g    = air gap length (m)
%  EE.ds   = slot depth (m)
%  EE.ws   = slot width (m)
%  EE.we   = end width (m)
%  EE.wc   = center width (m)
%  EE.wb   = base width (m)
%  EE.lc   = length of core (m)
% W        = structure of winding parameters
%  W.dw    = winding depth (m)
%  W.ww    = winding width (m)
% ksc      = scaling factor
% fn       = figure number
%
% Internal:
% wt       = total EE core width (m)
% ht       = total EE core height (m)
% xlbc     = position of lower left corner
% ylbc     = position of upper right corner
%
% Written by:
% G. Shane and S.D. Sudhoff
% School of Electrical and Computer Engineering
% 1285 Electrical Engineering Building
% West Lafayette, IN 47907-1285
% Phone: 765-494-3246
% Fax: 765-494-0676
% E-mail: sudhoff@ecn.purdue.edu
%
% Modified by Avyay Sah
% Email: asah@purdue.edu

% center coordinate
wt=2*EE.we+2*EE.ws+EE.wc;        % total width
ht=2*EE.wb+EE.g+2*EE.ds;       % total height
xlbc=-wt;
ylbc=-ht;
xhbc=xlbc;
yhbc=ylbc+2*(EE.wb+EE.ds);

% profile view
set(figure(fn),'color','white')
draw_EE_profileview(xlbc,ylbc,xhbc,yhbc,ksc,EE,W);
axis equal

end

% Profile view
function [] = draw_EE_profileview(xlbc,ylbc,xhbc,yhbc,ksc,EE,W)

% draw_EE_profileview  Draws a profile view of EE core inductor
%
% [] = draw_EE_profileview(xlbc,ylbc,ksc,EE,W)
%
% Inputs:
% xlbc     = x coordinate for left bottom corner (m)
% ykbc     = y coordinate for left bottom corner (m)
% ksc      = scaling factor
% EE       = structure of EI core parameters
%  EE.g    = air gap length (m)
%  EE.ds   = slot depth (m)
%  EE.ws   = slot width (m)
%  EE.we   = end width (m)
%  EE.wc   = center width (m)
%  EE.wb   = base width (m)
%  EE.lc   = length of core (m)
% W        = structure of winding parameters
%  W.dw    = winding depth (m)
%  W.ww    = winding width (m)
%
% Internal:
% xe*      = x-coordinates of e-core points (m)
% ye*      = y-coordinates of e-core points (m)
% xe       = vector of e-core point coordinates (m)
% ye       = vector of e-core point cordinates (m)
% xi*      = x-coordinates of i-core points (m)
% yi*      = y-coordinates of i-core points (m)
% xi       = vector of i-core point coordinates (m)
% yi       = vector of i-core point cordinates (m)
% xw*      = x-coordinates of winding points (m)
% yw*      = y-coordinates of winding points (m)
% xw       = vector of winding point coordinates (m)
% yw       = vector of winding point cordinates (m)
% dx       = displacement between widings (m)
%
% Written by:
% S.D. Sudhoff
% School of Electrical and Computer Engineering
% 1285 Electrical Engineering Building
% West Lafayette, IN 47907-1285
% Phone: 765-494-3246
% Fax: 765-494-0676
% E-mail: sudhoff@ecn.purdue.edu

% E-1 Core
xe0=xlbc;
xe1=xe0;
xe2=xe1+EE.we;
xe3=xe2;
xe4=xe3+EE.ws;
xe5=xe4;
xe6=xe5+EE.wc;
xe7=xe6;
xe8=xe7+EE.ws;
xe9=xe8;
xe10=xe9+EE.we;
xe11=xe10;

ye0=ylbc;
ye1=ye0+EE.ds+EE.wb;
ye2=ye1;
ye3=ye2-EE.ds;
ye4=ye3;
ye5=ye4+EE.ds;
ye6=ye5;
ye7=ye6-EE.ds;
ye8=ye7;
ye9=ye8+EE.ds;
ye10=ye9;
ye11=ye10-EE.ds-EE.wb;

xe=[xe0 xe1 xe2 xe3 xe4 xe5 xe6 xe7 xe8 xe9 xe10 xe11 xe0]*ksc;
ye=[ye0 ye1 ye2 ye3 ye4 ye5 ye6 ye7 ye8 ye9 ye10 ye11 ye0]*ksc;

fill(xe,ye,[0.75 0.75 0.75],'EdgeColor','k');
hold on

% E-2 Core
xee0=xhbc;
xee1=xee0;
xee2=xee1+EE.we;
xee3=xee2;
xee4=xee3+EE.ws;
xee5=xee4;
xee6=xee5+EE.wc;
xee7=xee6;
xee8=xee7+EE.ws;
xee9=xee8;
xee10=xee9+EE.we;
xee11=xee10;

yee0=yhbc;
yee1=yee0+EE.ds+EE.wb+2*EE.g;
yee2=yee1;
yee3=yee2-EE.ds;
yee4=yee3;
yee5=yee4+EE.ds-EE.g;
yee6=yee5;
yee7=yee6-EE.ds+EE.g;
yee8=yee7;
yee9=yee8+EE.ds;
yee10=yee9;
yee11=yee10-EE.ds-EE.wb-2*EE.g;

xee=[xee0 xee11 xee10 xee9 xee8 xee7 xee6 ...
    xee5 xee4 xee3 xee2 xee1 xee0]*ksc;
yee=[-yee0 -yee11 -yee10 -yee9 -yee8 -yee7...
    -yee6 -yee5 -yee4 -yee3 -yee2 -yee1 -yee0]*ksc;

fill(xee,yee,[0.75 0.75 0.75],'EdgeColor','k');

% winding
xw0=xe4-EE.csc;
xw1=xw0-W.ww;
xw2=xw1;

xw3=xw2+W.ww;
yw0=ye4+EE.cb;
yw1=yw0;
yw2=yw1+W.dw;
yw3=yw2;

xw=[xw0 xw1 xw2 xw3 xw0]*ksc;
yw=[yw0 yw1 yw2 yw3 yw0]*ksc;
 
fill(xw,yw,'y','EdgeColor','k');


xww0=xee4-EE.csc;
xww1=xww0-W.ww;
xww2=xww1;

xww3=xww2+W.ww;
yww0=yee4+EE.cb;
yww1=yww0;
yww2=yww1+W.dw;
yww3=yww2;

xww=[xww0 xww1 xww2 xww3 xww0]*ksc;
yww=[-yww0 -yww1 -yww2 -yww3 -yww0]*ksc;
 
fill(xww,yww,'y','EdgeColor','k');

xwww0=xe8-EE.cso;
xwww1=xwww0-W.ww;
xwww2=xwww1;

xwww3=xwww2+W.ww;
ywww0=ye8+EE.cb;
ywww1=ywww0;
ywww2=ywww1+W.dw;
ywww3=ywww2;

xwww=[xwww0 xwww1 xwww2 xwww3 xwww0]*ksc;
ywww=[ywww0 ywww1 ywww2 ywww3 ywww0]*ksc;
 
fill(xwww,ywww,'y','EdgeColor','k');

xwwww0=xee8-EE.cso;
xwwww1=xwwww0-W.ww;
xwwww2=xwwww1;

xwwww3=xwwww2+W.ww;
ywwww0=yee8+EE.cb;
ywwww1=ywwww0;
ywwww2=ywwww1+W.dw;
ywwww3=ywwww2;

xwwww=[xwwww0 xwwww1 xwwww2 xwwww3 xwwww0]*ksc;
ywwww=[-ywwww0 -ywwww1 -ywwww2 -ywwww3 -ywwww0]*ksc;
 
fill(xwwww,ywwww,'y','EdgeColor','k');

%plot(xw,yw)
fill(xw,yw,'y','EdgeColor','k');
fill(xww,yww,'y','EdgeColor','k');
fill(xwww,ywww,'y','EdgeColor','k');
fill(xwwww,ywwww,'y','EdgeColor','k');
hold off

end