function pmac_plotting_v0_2(D,G,W,R,S,M,C,fn)
% pmac_plotting_v0_2 is plotting function that create a 2D axial plot,
%                    2D radial plot, and 3D plot of a pmac machine.
%
% Call:
% pmac_plotting_v0_2(D,G,W,R,S,M,C,fn)
%
% Inputs:
% D        =   design specificatin tructure
% G        =   geometry parameter structure
% W        =   winding parameter structure
% R        =   rotor material parameter structure
% S        =   stator material parameter structure
% M        =   magnet material parameter structure
% C        =   conductor material parameter structure
% fn       =   1x3 figure number vector i.e [1 2 3]
%
% Date/Version:
% v0.1     =   first draft version.
% v0.2     =   corrected axial l2 axis variable and alphatt bug.
%              Added slot feed to end windings in axial plot.
% Notes:
% Internal function revisions are reported in this main function header
%
% Written by:
% Jacob Krizan
% Purdue University - September 2011
% jkrizan@purdue.edu

% Append parameters
G.alphas = 1-G.alphatt;
G.alphai = 0;
G.di = 0;

% Radial Machine Plot
figure(fn(1))
clf
PMSM_radial_plot(D,G,W,R,S,M,C);
axis square

% Axial Machine Plot
figure(fn(2))
clf
PMSM_axial_plot(D,G,W,R,S,M,C);
axis square

% 3D Machine Plot
figure(fn(3))
clf
axis off
hold on
PMSM_3D_plot(D,G,R,S,M,C);
%End pmac_plotting_v1_0----------------------------------------------------

function PMSM_radial_plot(D,G,W,R,S,M,C)
%  FUNCTION:    PMSM_radial_plot(D,G,W,R,S,M,C);
%
%  DESCRIPTION: function to plot radial view of PMSM
%
%  INPUTS:      E.x         - Electrical Parameters
%               G.x         - Geometry Parameters
%               M.x         - Material Parameters
%               
%  OUTPUTS:     Depiction of 2D PMSM radial view on open figure
%
%  AUTHOR:      Jacob Krizan
%               Purdue University - June 2011
%               jkrizan@purdue.edu
%

%% Material Colors
scolor = colormap([0.0742    0.3633    0.5586]);
sdesc = S.desc;
rcolor = colormap([0.4336    0.5117    0.6289]);
rdesc = R.desc;
pmcolor = colormap([0.1992    0.2422    0.3086]);
pmdesc = M.desc;
if C.row == 8890
    ccolor = colormap([0.9531    0.3008    0.1719]);
elseif C.row == 2705
    ccolor = colormap([0.7305    0.7578    0.8008]);
else
    ccolor = colormap([0.9531    0.3008    0.1719]);
end
cdesc = C.desc;

%% Plot colors for legend sake
hold on
draw_box_2D(1e-3,1e-3,0,0,0,scolor);
draw_box_2D(1e-3,1e-3,0,0,0,rcolor);
draw_box_2D(1e-3,1e-3,0,0,0,ccolor);
draw_box_2D(1e-3,1e-3,0,0,0,pmcolor);

%% Plot casing
ric = G.rri+G.drb+G.dm+G.g+G.dst+G.dsb;
dc = .02*ric;
draw_tube_2D(ric,ric+dc,0,100,'k',0,0);

%% Plot Backiron
ri = G.rri+G.drb+G.dm+G.g+G.dst;
draw_tube_2D(ri,ri+G.dsb,0,100,scolor,0,0);

%% Plot Stator Teeth Box
ri = G.rri+G.drb+G.dm+G.g;
tooth_percent = (1-G.alphas)*100/G.Ss;
h =2*pi*ri*tooth_percent/100;
w = G.dst;
xoff = ri+G.dst/2;
yoff = 0;
section_angle = 360/G.Ss; % (deg)
for i = 1:G.Ss
    a = (i-1)*section_angle-section_angle/2; % (deg)
    draw_box_2D(w,h,a,xoff,yoff,scolor);
end

%% Plot Conductors

% Define Radii and Constants
r = G.rri+G.drb+G.dm+G.g;
rc = sqrt(W.ac/pi);
rl = r+G.dst-rc;
ri = G.rri+G.drb+G.dm+G.g;
tooth_percent = (1-G.alphas)*100/G.Ss;
tooth_width =2*pi*ri*tooth_percent/100;
fraction = 1/G.Ss;

% Loop for each slot
section_angle = 360/G.Ss; % (deg)
for j = 1:G.Ss
    % Define Angle
    a = (j-1)*section_angle; % (deg)

    % Initialize
    Nt = W.Ns(j);
    e = zeros(1,Nt);
    xoff = zeros(Nt,1);
    yoff = zeros(Nt,1);
    xofft = zeros(Nt,1);
    yofft = zeros(Nt,1);
    
    % Compute Layer Widths
    Lw = 0;
    Nwire = 0;
    rlw = 0;
    i = 1;
    while Nt > 0
        rlw(i) = rl-2*(i-1)*rc;
        Lw(i) = 2*pi*(rlw(i))*fraction-tooth_width;
        Nwire(i) = floor(Lw(i)/(2*rc));
        % Only use remaining turns in layer
        if Nwire(i) > Nt
            Nwire(i) = Nt;
        end
        % Track total wires
        if Nwire(i) <= 0 && Nt ~= 0 
            % Deal with wires that do not fit
            Nwire(i) = 1;
            Nt = Nt-1;
            e(length(e)-Nt) = 1;
        else
            % Subtract Layer wires from Nt
            Nt = Nt-Nwire(i);
        end
        i = i+1;    % Loop counter
    end

    % Place Wires in Layers (find offset positions)
    Nlayers = length(Lw);
    idx = 0;
    for i=1:Nlayers
        sp = Lw(i)/Nwire(i);
        for j = 1:Nwire(i)
            xoff(j+idx) = rl-2*(i-1)*rc;
            yoff(j+idx) = j*sp-Lw(i)/2-sp/2;
        end
        idx = idx +Nwire(i);
    end
    
    % Rotate Offsets
    Q = [cosd(a) -sind(a);...
         sind(a) cosd(a)];
    for i = 1:length(xoff)
        xy = Q*[xoff(i) ; yoff(i)];
        xofft(i) = xy(1);
        yofft(i) = xy(2);
    end

    % Plot Conductors
    for i = 1:length(xofft)
        % Choose Color
        if e(i) == 1
            ccolor = 'r';
        end
        % Draw single conductor
        draw_tube_2D(0,rc,0,100,ccolor,xofft(i),yofft(i),20);
        draw_tube_2D(.8*rc,rc,0,100,'k',xofft(i),yofft(i),20);
    end
end

%% Plot Q-axis Rotor Teeth
ri = G.rri+G.drb;
tooth_percent = G.alphai*100/G.P;
pole_angle = 360/G.P; % (deg)
for i = 1:G.Ss
    tooth_angle = (i-1)*pole_angle-G.alphai*pole_angle/2; % (deg)
    draw_tube_2D(ri,ri+G.di,tooth_angle,tooth_percent,rcolor,0,0);  
end

%% Plot Magnets

% 1st half
ri = G.rri+G.drb;
pole_angle = 360/G.P; % (deg)
mag_percent = (G.alphapm)*100/(2*G.P);
for i = 1:G.Ss
    mag_angle = (i-1)*pole_angle+(1-G.alphapm)*pole_angle/2;
    draw_tube_2D(ri,ri+G.dm,mag_angle,mag_percent,pmcolor,0,0);  
end

% 2nd half
ri = G.rri+G.drb;
pole_angle = 360/G.P; % (deg)
mag_percent = (G.alphapm)*100/(2*G.P);
for i = 1:G.Ss
    mag_angle = (i-1)*pole_angle+pole_angle/2;
    draw_tube_2D(ri,ri+G.dm,mag_angle,mag_percent,pmcolor,0,0);  
end

%% Plot rotor backiron
ri = G.rri;
draw_tube_2D(ri,ri+G.drb,0,100,rcolor,0,0);

%% Plot Spacer
c = colormap([0.992 0.930 0.836]);
draw_tube_2D(D.rrs,G.rri,0,100,c,0,0);

%% Plot Shaft
draw_tube_2D(0,D.rrs,0,100,'k',0,0);

%% Legend, Title and the like

ylabel('m')
title('PMSM Radial Cross-Section')
legend(sdesc,rdesc,cdesc,pmdesc)
hold off
axis([-1.05*(ric+dc) 1.05*(ric+dc) -1.05*(ric+dc) 1.05*(ric+dc)])
%End PMSM_radial_plot------------------------------------------------------

function PMSM_axial_plot(D,G,W,R,S,M,C)
%  FUNCTION:    PMSM_axial_plot(D,G,W,R,S,M,C);
%
%  DESCRIPTION: function to plot axial view of PMSM
%
%  INPUTS:      E.x         - Electrical Parameters
%               G.x         - Geometry Parameters
%               M.x         - Material Parameters
%               
%  OUTPUTS:     Depiction of PMSM on open figure
%
%  AUTHOR:      Jacob Krizan
%               Purdue University - June 2011
%               jkrizan@purdue.edu

%% Material Colors
scolor = colormap([0.0742    0.3633    0.5586]);
sdesc = S.desc;
rcolor = colormap([0.4336    0.5117    0.6289]);
rdesc = R.desc;
pmcolor = colormap([0.1992    0.2422    0.3086]);
pmdesc = M.desc;
if C.row == 8890
    ccolor = colormap([0.9531    0.3008    0.1719]);
elseif C.row == 2705
    ccolor = colormap([0.7305    0.7578    0.8008]);
else
    ccolor = colormap([0.9531    0.3008    0.1719]);
end
cdesc = C.desc;

%% Plot colors for legend sake
hold on
draw_box_2D(1e-3,1e-3,0,0,0,scolor);
draw_box_2D(1e-3,1e-3,0,0,0,rcolor);
draw_box_2D(1e-3,1e-3,0,0,0,ccolor);
draw_box_2D(1e-3,1e-3,0,0,0,pmcolor);

%% Plot casing
ric = G.rri+G.drb+G.dm+G.g+G.dst+G.dsb;
dc = .02*ric;      % depth of casing
draw_box_2D(2*(ric+dc),G.l+2*(W.lew+D.leo+dc),0,0,0,'k');
draw_box_2D(2*(ric),G.l+2*W.lew+2*D.leo,0,0,0,'w');

%% Plot Stator
ri = G.rri+G.drb+G.dm+G.g;
draw_box_2D(G.dsb+G.dst,G.l,0,ri+(G.dsb+G.dst)/2,0,scolor);
draw_box_2D(G.dsb+G.dst,G.l,0,-ri-(G.dsb+G.dst)/2,0,scolor);

%% Plot End Windings
ri = G.rri+G.drb+G.dm+G.g;
dl = 0.025*W.lew;   % depth of end winding outline

draw_box_2D(G.dst,D.leo,0,ri+G.dst/2,(G.l+D.leo)/2,ccolor);
draw_box_2D(G.dst,D.leo,0,-ri-G.dst/2,(G.l+D.leo)/2,ccolor);
draw_box_2D(G.dst,D.leo,0,ri+G.dst/2,-(G.l+D.leo)/2,ccolor);
draw_box_2D(G.dst,D.leo,0,-ri-G.dst/2,-(G.l+D.leo)/2,ccolor);

draw_box_2D(G.dst,W.lew,0,ri+G.dst/2,(G.l+W.lew)/2+D.leo,'k');
draw_box_2D(G.dst,W.lew,0,-ri-G.dst/2,(G.l+W.lew)/2+D.leo,'k');
draw_box_2D(G.dst,W.lew,0,-ri-G.dst/2,-(G.l+W.lew)/2-D.leo,'k');
draw_box_2D(G.dst,W.lew,0,ri+G.dst/2,-(G.l+W.lew)/2-D.leo,'k');

draw_box_2D(G.dst-dl,W.lew-dl,0,ri+G.dst/2,(G.l+W.lew)/2+D.leo,ccolor);
draw_box_2D(G.dst-dl,W.lew-dl,0,-ri-G.dst/2,(G.l+W.lew)/2+D.leo,ccolor);
draw_box_2D(G.dst-dl,W.lew-dl,0,-ri-G.dst/2,-(G.l+W.lew)/2-D.leo,ccolor);
draw_box_2D(G.dst-dl,W.lew-dl,0,ri+G.dst/2,-(G.l+W.lew)/2-D.leo,ccolor);


%% Plot Magnets
ri = G.rri+G.drb;
draw_box_2D(G.dm,G.l,0,ri+G.dm/2,0,pmcolor);
draw_box_2D(G.dm,G.l,0,-ri-G.dm/2,0,pmcolor);

%% Plot rotor backiron
ri = G.rri;
draw_box_2D(G.drb,G.l,0,ri+G.drb/2,0,rcolor);
draw_box_2D(G.drb,G.l,0,-ri-G.drb/2,0,rcolor);

%% Plot Spacer
c = colormap([0.992 0.930 0.836]);
draw_box_2D(2*G.rri,G.l,0,0,0,c);

%% Plot Shaft
lsh = G.l+D.lfs+D.lbs;
draw_box_2D(2*D.rrs,lsh,0,0,0,'k');

%% Legend, Titles, and the like
ylabel('m')
title('PMSM Axial Cross-Section')
legend(sdesc,rdesc,cdesc,pmdesc)
hold off

% Ensure proper axis
l1 = G.l/2+(W.lew+D.leo+dc);
l2 = ric+dc;
if l1 > l2
    axis([-1.05*l1 1.05*l1 -1.05*l1 1.05*l1])
else
    axis([-1.05*l2 1.05*l2 -1.05*l2 1.05*l2])
end
%End PMSM_axial_plot-------------------------------------------------------

function PMSM_3D_plot(D,G,R,S,M,C)
%  FUNCTION:    PMSM_3D_plot(D,G,R,S,M,C);
%
%  DESCRIPTION: function to plot 3D view of PMSM
%
%  INPUTS:      E.x         - Electrical Parameters
%               G.x         - Geometry Parameters
%               M.x         - Material Parameters
%               
%  OUTPUTS:     Depiction of PMSM on open figure
%
%  AUTHOR:      Jacob Krizan
%               Purdue University - June 2011
%               jkrizan@purdue.edu


%% Material Colors
scolor = colormap([0.0742    0.3633    0.5586]);
sdesc = S.desc;
rcolor = colormap([0.4336    0.5117    0.6289]);
rdesc = R.desc;
pmcolor = colormap([0.1992    0.2422    0.3086]);
pmdesc = M.desc;
if C.row == 8890
    ccolor = colormap([0.9531    0.3008    0.1719]);
elseif C.row == 2705
    ccolor = colormap([0.7305    0.7578    0.8008]);
else
    ccolor = colormap([0.9531    0.3008    0.1719]);
end
cdesc = C.desc;

%% Plot casing
ric = G.rri+G.drb+G.dm+G.g+G.dst+G.dsb;
dc = .02*ric;
c = colormap([0.15 0.15 0.15]);
draw_tube_3D(ric,ric+dc,0,100,G.l,c,'-',1,0,0,0,0,100);

%% Plot Backiron
ri = G.rri+G.drb+G.dm+G.g+G.dst;
draw_tube_3D(ri,ri+G.dsb,0,100,G.l,scolor,'-',1,0,0,0,0,100);

%% Plot Stator Teeth Box
ri = G.rri+G.drb+G.dm+G.g;
tooth_percent = (1-G.alphas)*100/G.Ss;
h =2*pi*ri*tooth_percent/100;
w = G.dst;
yoff = 0;
section_angle = 360/G.Ss; % (deg)
for i = 1:G.Ss
    a = (i-1)*section_angle-section_angle/2; % (deg)
    xoff = (ri+w/2)*cosd(a);
    zoff = (ri+w/2)*sind(a);
    draw_box_3D(G.l,w,h,xoff,yoff,zoff,a,2,[1 1 1 1 1 1],1,'-',scolor);
end

%% Plot rotor backiron
ri = G.rri;
draw_tube_3D(ri,ri+G.drb,0,100,G.l,rcolor,'-',1,0,0,0,0,100);

%% Plot Magnets
ri = G.rri+G.drb;
pole_angle = 360/G.P; % (deg)
mag_percent = 2*(G.alphapm)*100/(2*G.P);
for i = 1:G.Ss
    mag_angle = (i-1)*pole_angle+(1-G.alphapm)*pole_angle/2;
    draw_tube_3D(ri,ri+G.dm,mag_angle,mag_percent,G.l,pmcolor,'-',1,0,0,0,0,10);
end

%% Plot Rotor Teeth
if G.alphai ~= 0
    ri = G.rri+G.drb;
    tooth_percent = G.alphai*100/G.P;
    pole_angle = 360/G.P; % (deg)
    for i = 1:G.Ss
        tooth_angle = (i-1)*pole_angle-G.alphai*pole_angle/2; % (deg)
        draw_tube_3D(ri,ri+G.di,tooth_angle,tooth_percent,G.l,rcolor,'-',1,0,0,0,0,10);
    end
end
%% Plot Spacer
c = colormap([0.992 0.930 0.836]);
draw_tube_3D(D.rrs,G.rri,0,100,G.l,c,'-',1,0,0,0,0,100);

%% Plot Shaft
c = colormap([0.15 0.15 0.15]);
draw_tube_3D(0,D.rrs,0,100,G.l+D.lfs+D.lbs,c,'-',1,0,0,0-D.lfs,0,100);

% Set figure properties
view(3)
axis equal
%End PMSM_3D_plot----------------------------------------------------------

function draw_box_2D(w,h,a,xoff,yoff,c)
% draw_box draws a 2D box as specified by the inputs
%
% draw_box_2D(w,h,a,xoff,yoff,c);
%
% NOTE: Must open figure and put hold on before calling 
%
% Inputs:
%   w       = width of box
%   h       = height of box
%   a       = angle of box for rotation around 0,0 (degrees)
%   xoff    = x axis offset of box
%   yoff    = y axis offset of box
%   c       = string denoting color (i.e. 'r')
%
% Outputs:
%   plot of 2D tube drawing on open figure
%
% Written by:
% Jacob Krizan
% Purdue University - 2010
% E-mail: jkrizan@purdue.edu

% Vertices of rectangular box
v1 = [-w/2 + xoff , -h/2 + yoff ];
v2 = [w/2  + xoff , -h/2 + yoff ];
v3 = [w/2  + xoff ,  h/2 + yoff ];
v4 = [-w/2 + xoff ,  h/2 + yoff ];

% Create Rotation Matrix
Q = [cosd(a) -sind(a);...
     sind(a) cosd(a)];

% Rotate Verticies
v1 = Q*v1';
v2 = Q*v2';
v3 = Q*v3';
v4 = Q*v4';

% Compile Vectors and add offset
x = [v1(1) v2(1) v3(1) v4(1)];
y = [v1(2) v2(2) v3(2) v4(2)];

% Plot Drawing
fill(x,y,c,'EdgeAlpha',0)
%End draw_box_2D-----------------------------------------------------------

function draw_box_3D(w,d,h,xoff,yoff,zoff,theta,raxis,fv,tp,ls,c)
% draw_box draws a 3D box as specified by the inputs
%
% draw_box_3D(w,d,h,xoff,yoff,zoff,theta,raxis,fv,tp,ls,c);
%
% NOTE: Must open figure and put hold on before calling 
%
% NOTE: (0,0,0) is in the center of the box
%
% Inputs:
%   w       = width of box
%   d       = depth of box
%   h       = height of box
%   xoff    = x axis offset of box
%   yoff    = y axis offset of box
%   zoff    = z axis offset of box
%   theta   = rotation angle around center point of box (degrees)
%   raxis   = rotation axis (1=x , 2=y , 3=z)
%   fv      = vector of faces to draw (i.e. [1 1 1 1 1 1]) 
%   tp      = transparency of box faces (0 to 1)
%   ls              = string denoting line style (i.e. ':')
%   c       = string denoting color (i.e. 'r')
%
% Outputs:
%   plot of 3D tube drawing on open figure
%
% Written by:
% Jacob Krizan
% Purdue University - June 2011
% E-mail: jkrizan@purdue.edu

%% Vertices of rectangular box
v1 = [-d/2,-w/2,-h/2];
v2 = [ d/2,-w/2,-h/2];
v3 = [ d/2, w/2,-h/2];
v4 = [-d/2, w/2,-h/2];
v5 = [-d/2,-w/2, h/2];
v6 = [ d/2,-w/2, h/2];
v7 = [ d/2, w/2, h/2];
v8 = [-d/2, w/2, h/2];

%% Add Rotation

% Rotation about x-axis
if raxis == 1
    rotate = [       1           0            0
                     0      cosd(theta) -sind(theta)
                     0      sind(theta)  cosd(theta)];
% Rotation about y-axis
elseif raxis == 2
    rotate = [ cosd(theta)       0       sind(theta)
                     0           1            0
              -sind(theta)       0       cosd(theta)];
% Rotation about z-axis
elseif raxis == 3
    rotate = [ cosd(theta) -sind(theta)       0
               sind(theta)  cosd(theta)       0
                    0            0            1     ];
% Error reporting
else
    disp('Error: Rotational axis parameter needs to be 1, 2, or 3')
    return
end
    
% Rotate vertices
v1 = v1*rotate;
v2 = v2*rotate;
v3 = v3*rotate;
v4 = v4*rotate;
v5 = v5*rotate;
v6 = v6*rotate;
v7 = v7*rotate;
v8 = v8*rotate;

%% Add offset
vx = [v1(1); v2(1); v3(1); v4(1); v5(1); v6(1); v7(1); v8(1)] + xoff;
vy = [v1(2); v2(2); v3(2); v4(2); v5(2); v6(2); v7(2); v8(2)] + yoff; 
vz = [v1(3); v2(3); v3(3); v4(3); v5(3); v6(3); v7(3); v8(3)] + zoff; 

%% Compile vertices
%vert = [v1; v2; v3; v4; v5; v6; v7; v8];
vert = [vx vy vz];


%% Faces of rectangular box
f1 = [1 2 3 4];
f2 = [3 4 8 7];
f3 = [4 1 5 8];
f4 = [1 2 6 5];
f5 = [2 3 7 6];
f6 = [5 6 7 8];
faces = [f1; f2; f3; f4; f5; f6];

% Compile Face List
for n = 1:6
    if fv(n) == 0
        faces(n,:) = [1 1 1 1];
    end
end

%% Plot Drawing
patch('vertices',vert,'faces',faces,'facecolor',c,...
    'FaceAlpha',tp,'LineStyle',ls)
%End draw_box_3D-----------------------------------------------------------

function draw_tube_2D(ri,ro,start_degree,percent,color,xoff,yoff,pts)
% draw_tube_2D draws a 2D tube as specified by the inputs
%
% draw_tube_2D(ri,ro,start_degree,percent,color,xoff,yoff);
%
% NOTE: Must open figure and put hold on before calling 
%
% Inputs:
%   ri              = inner radius of tube
%   ro              = outer radius of tube
%   start_degree    = degree to start tube (0 to 360) (degrees)
%   percent         = percent of tube to draw (0 to 100) (percent)
%   color           = string denoting color (i.e. 'r')
%   xoff            = x axis offset of tube
%   yoff            = y axis offset of tube
%
% Outputs:
%   plot of 2D tube drawing on open figure
%
% Written by:
% Jacob Krizan
% Purdue University - 2009
% E-mail: jkrizan@purdue.edu
%

% Assign number of points
if nargin == 8
    points = pts;
else
    points = 5e3*(percent/100);
end

% Construct Tube Edges
off = start_degree*(pi/180);        % offset to start circle
theta1 = linspace(off,(percent/100)*2*pi+off,points); 
x1 = ro*cos(theta1);
y1 = ro*sin(theta1);   
theta2 = linspace((percent/100)*2*pi+off,off,points); 
x2 = ri*cos(theta2);
y2 = ri*sin(theta2);           

% Compile Vector of edge points
x = [x1 x2];
y = [y1 y2];

% Plot tube on open figure
fill(x+xoff,y+yoff,color,'EdgeAlpha',0)
hold on
%End draw_tube_2D----------------------------------------------------------

function draw_tube_3D(ri,ro,start_degree,percent,l,color,ls,transp,rotate,...
    xoff,yoff,zoff,points)
% draw_tube draws a 3D tube as specified by the inputs
%
% draw_tube_3D(ri,ro,start_degree,percent,l,color,ls,transp,rotate,...
%              xoff,yoff,zoff,points);
%
% NOTE: Must open figure and put hold on before calling 
%
% Inputs:
%   ri              = inner radius of tube
%   ro              = outer radius of tube
%   start_degree    = degree to start tube (0 to 360) (degrees)
%   percent         = percent of tube to draw (0 to 100) (percent)
%   l               = length of tube
%   color           = string denoting color (i.e. 'r')
%   ls              = string denoting line style (i.e. ':')
%   transp          = transparency of tube face (0 to 1)
%   rotate          = set 1 to rotate by 90 degrees to 0 for no rotation
%   xoff            = x axis offset of tube
%   yoff            = y axis offset of tube
%   zoff            = z axis offset of tube
%   points          = number of dicrete points that discribe circle sector
%
% Outputs:
%   plot of 3D tube drawing on open figure
%
% Written by:
% Jacob Krizan
% Purdue University - June 2011
% E-mail: jkrizan@purdue.edu

%% Check for errors
if percent > 100
    disp('Error: Percent cannot be great than 100%!')
    return
end

if points < 4
    disp('Error: Must use minimum of four points to discribe sector!')
    return
end

%% Draw Outer and Inner Faces of Tubes-------------------------------------
r = [ri ro];                            % define radius vector
% draw first inner, then outer tube face
for n = 1:2
    % construct Tube
    off = start_degree*(pi/180);        % offset to start circle
    theta = linspace(off,(percent/100)*2*pi+off,points); 
    x = r(n)*cos(theta);
    y = linspace(-l/2,l/2,points);
    zv = r(n)*sin(theta);               % z height vector
    z = meshgrid(zv,ones(length(y),1)); % make z into matrix for surf 
    
    % rotate tube if specified by input
    if rotate ~= 0
        xp = x;
        x = y;
        y = xp;
        z = z'; 
    end
    
    % plot outer and inner faces of tube
    surf(x+xoff,y+yoff,z+zoff,'FaceColor',color,'FaceAlpha',...
        transp,'EdgeAlpha',0)
end
 
%% Draw End Faces of Tubes-------------------------------------------------
% get vector of y edge points
x = [r(1)*cos(theta) r(2)*cos(theta)]+xoff;
z = [r(1)*sin(theta) r(2)*sin(theta)]+zoff;
L = length(x);

% build vertices and face arrays
vert1 = zeros(L,3);
vert2 = zeros(L,3);
face = zeros(1,L);
for i = 1:L
    % rotate tube if specified by input
    if rotate == 0
        % compile verticies of end face
        vert1(i,:) = [x(i) -l/2+yoff z(i)];
        vert2(i,:) = [x(i) l/2+yoff z(i)];
    else
        % compile verticies of end face
        vert1(i,:) = [-l/2+yoff x(i) z(i)];
        vert2(i,:) = [l/2+yoff x(i) z(i)];    
    end
    if i > L/2
        % compile end face out of outer verticies
        face(i) = 1.5*L+1-i;
    else
        % compile end face out of inner verticies
        face(i) = i; 
    end
end

% plot end faces of tube
patch('vertices',vert1,'faces',face,'facecolor',color,...
    'FaceAlpha',transp,'LineStyle',ls)
patch('vertices',vert2,'faces',face,'facecolor',color,...
    'FaceAlpha',transp,'LineStyle',ls)

%% Draw Tube Sector Cut Ends
if percent < 100
    %% Draw first cut end
    
    % Generate x and z coordinates     
    x = [r(1)*cos(off) r(2)*cos(off)]+xoff;
    z = [r(1)*sin(off) r(2)*sin(off)]+zoff;
        
    % Constuct Corners
    if rotate == 0
        x1 = [x(1),yoff - l/2,z(1)];
        x2 = [x(2),yoff - l/2,z(2)];
        x3 = [x(2),yoff + l/2,z(2)];
        x4 = [x(1),yoff + l/2,z(1)];
    else
        x1 = [yoff - l/2,x(1),z(1)];
        x2 = [yoff - l/2,x(2),z(2)];
        x3 = [yoff + l/2,x(2),z(2)];
        x4 = [yoff + l/2,x(1),z(1)];
    end
    
    % Compile corners and face
    corners1 = [x1 ; x2 ; x3 ; x4];
    face1 = [1 2 3 4];
        
    % Plot first cut end
    patch('vertices',corners1,'faces',face1,'facecolor',color,...
    'FaceAlpha',transp,'LineStyle',ls)
    
    %% Draw second cut end
    
    % Generate x and z coordinates     
    x = [r(1)*cos(off + (percent/100)*2*pi) r(2)*cos(off + (percent/100)*2*pi)]+xoff;
    z = [r(1)*sin(off + (percent/100)*2*pi) r(2)*sin(off + (percent/100)*2*pi)]+zoff;
        
    % Constuct Corners
    if rotate == 0
        x1 = [x(1),yoff - l/2,z(1)];
        x2 = [x(2),yoff - l/2,z(2)];
        x3 = [x(2),yoff + l/2,z(2)];
        x4 = [x(1),yoff + l/2,z(1)];
    else
        x1 = [yoff - l/2,x(1),z(1)];
        x2 = [yoff - l/2,x(2),z(2)];
        x3 = [yoff + l/2,x(2),z(2)];
        x4 = [yoff + l/2,x(1),z(1)];
    end
    
    % Compile corners and face
    corners1 = [x1 ; x2 ; x3 ; x4];
    face1 = [1 2 3 4];
        
    % Plot first cut end
    patch('vertices',corners1,'faces',face1,'facecolor',color,...
    'FaceAlpha',transp,'LineStyle',ls)

end
%End draw_tube_3D----------------------------------------------------------
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                