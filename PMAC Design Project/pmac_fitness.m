function f=pmac_fitness(p,D,fn)
% pmac_fitness is a fitness function for a design of a surface
%              mounted permanent magnet ac machine as used in Chapter 9,
%              Section 11 of "Power Magnetic Devices: A Multi-Objective
%              Design Approach," by S.D. Sudhoff
% 
% Call:
% f=pmac_fitness(p,D)
% f=pmac_fitness(p,D,fn)
% 
% Inputs:
% p        = vector of design variables
%  p(1)    = stator steel code (an integer)
%  p(2)    = rotor steel code (an integer)
%  p(3)    = conductor code (an integer)
%  p(4)    = magnet code (an integer)
%  p(5)    = number of pole pairs
%  p(6)    = depth of magnetically inert region (m)
%  p(7)    = depth of rotor back iron (m)
%  p(8)    = depth of magnet (m)
%  p(9)    = airgap (m)
%  p(10)   = depth of tooth base (m)
%  p(11)   = tooth tip fraction
%  p(12)   = depth of stator backiron (m)
%  p(13)   = permanent magnet fraction
%  p(14)   = active length (m)
%  p(15)   = pk phs cond.density(cond/rad)
%  p(16)   = relative 3rd harmon cond. density
%  p(17)   = Kt1: torque multiplier
%  p(18)   = Kt2: torque multiplier
%  p(19)   = Kt3: torque multiplier
%  p(20)   = Kt4: torque multiplier
%  p(21)   = Kt5: torque multiplier
%  p(22)   = Kt6: torque multiplier
% D        = vector of design specifications
%  D.vdc   = dc voltage (V)
%  D.wrm   = mechanical rotor speed (rad/s) 
%  D.Te    = req. torque in motor mode(Nm)
%  D.rrs   = shaft radius
%  D.kpf   = packing factor                                     
%  D.leo   = end winding offset (W)
%  D.vfs   = forward semiconductor drop (V)
%  D.J     = number of rotor positions
%  D.mlim  = mass limit (kg)
%  D.Pllim = loss limit (W)
%  D.atar  = maximum tooth aspect ratio
%  D.aso   = slot opening factor
%  D.nspp  = slots per pole per phase
%  D.km    = multiple of magnet Hci to get Hlim
%  D.lfs   = length of front shaft (m) 
%  D.lbs   = length of back shaft (m) 
% fn       = figure number (and flag to report)
%
% Outputs:
% sd number of signifcant digits
% f fitness vector (if fn is not specified)
%  f(1)    = reciprical of mass (1/m^3)
%  f(2)    = reciprical of loss (1/kg)
% f structure of selected design variables (if fn=0)
%  f.m     = mass (kg)
%  f.mss   = mass of sttor steel (kg)
%  f.mrs   = mass of rotor steel (kg)
%  f.mpm   = mass of permanent magnet (kg)
%  f.mcd   = mass of conductor (kg)
%  f.Pr    = resistive loss (W)
%  f.Ps    = semiconductor conduction loss (W)
%  f.Pc    = core loss (W)
%  f.Pl    = total loss (W)
%  f.P     = number of poles
%  f.Is    = rms current (A)
%  f.J     = rms current density (A/m^2)
%  f.ac    = conductor cross sectional area (m^2)
%  f.st    = stator steel code
%  f.rt    = rotor steel code
%  f.ct    = conductor code
%  f.mt    = magnet code
%  f.aslt  = slot area (m^3)
%  f.l     = active length (m)
%  f.rss   = radius of stator shell (m)
%  f.alphat= tooth fraction
%  f.alphapm=permanent magnet fraction
% 
% Internal:
% sixrt2opi= handy constant (6sqrt(2)/pi)
% NC       = number of constraints
% CS       = number of constraints satisfied
% CI       = number of constraints imposed
% c*       = constraint variable (1 if satisfied, 0 if not)
% S        = structure of stator steel parameters
%  S.desc  = steel description
%  S.Blim  = recommended limit on B to avoid saturation (T)
%  S.row   = mass density (kg/m^3)
%  S.k     = thermal conductivity (W/(m*K))
%  S.c     = specific heat capacity (J/(kg*K))
% R        = structure of rotor steel parameters
%  R.desc  = steel description
%  R.Blim  = recommended limit on B to avoid saturation (T)
%  R.row   = mass density (kg/m^3)
%  R.k     = thermal conductivity (W/(m*K))
%  R.c     = specific heat capacity (J/(kg*K))
% C        = structure of conductor parameters (see pmac_conductor_catalog)
% M        = structure of magnet parameters
%  M.desc  = magnet desription
%  M.br    = residual flux density (T)
%  M.x     = susceptibility 
%  M.Hlim  = minimum field intensity before demag (A/m)
%           (Taken to be 1/2 the intrinsic coercive force at room temp.)
%  M.row   = mass density (kg/m^3)
% G        = structure with machine geometry (see pmac_geometry)
% W        = structure of winding parameters (see pmac_winding)
% I        = structure of current (see pmac_current)
% mss      = mass of stator steel (kg)
% mrs      = mass of rotor steel (kg)
% mcd      = mass of conductor (kg)
% m        = total mass (kg)
% E        = structure of lump. parameters (see pmac_electrical_parameters)
% wr       = electrical rotor speed (rad/s)
% vqsr     = q-axis voltage (V)
% vdsr     = d-axis voltage (V)
% vpkll    = peak line-to-line voltage (V)
% Te       = torque (Nm)
% I0       = structure of current excitation data for zero 
%            current conditions (see pmac_current)
% F0       = structure of field quantities for zero current conditions 
%            (see pmac_field_analysis)
% F        = structure of field quantities for normal conditions
%            (see pmac_field_analysis)
% Pr       = resistive loss in machine (W)
% Ps       = semiconductor conduction loss (W)
% Pl       = total loss (W)
% Tec      = corrected torque (corrected for hysteresis loss) (Nm)
% Pc       = core loss (W)
% Plagg    = aggregate power loss (W)
% Pinaggtot= Total Aggregate Input (W)
% Poutagg  = Aggregate Output (W)
% efficiency= Aggregate Efficiency
%
% Written by:
% S.D. Sudhoff
% School of Electrical and Computer Engineering
% 1285 Electrical Engineering Building
% West Lafayette, IN 47907-1285
% Phone: 765-494-3246
% Fax: 765-494-0676
% E-mail: sudhoff@ecn.purdue.edu
%
% Modified by Avyay Sah
% Email: asah@purdue.edu

% handy number
persistent sixrt2opi;
if isempty(sixrt2opi)
    sixrt2opi=6*sqrt(2)/pi;
end

% Constraints
NC=18;                             % number of constraints
  
% Material Types-----------------------------------------------------------
S=steel_catalog(p(1));             % stator steel
R=steel_catalog(p(2));             % rotor steel
C=pmac_conductor_catalog(p(3));    % conductor
M=pmac_magnet_catalog(p(4));       % magnet

% Geometry-----------------------------------------------------------------
G.P=round(p(5))*2;
G.rrs=D.rrs;                         
G.di=p(6);   
G.drb=p(7);
G.dm=p(8);  
G.g=p(9);
G.dtb=p(10);
G.dttc=0;
G.dtte=0;
G.alphat=p(11);
G.alphatt=G.alphat;
G.dsb=p(12);
G.alphapm=p(13);
G.l=p(14);
G.Ss=3*G.P*D.nspp;
G.phiss1=pi/G.Ss;
G=pmac_geometry(G);

c1=lte(G.dst/G.wtb,D.atar); 

% Winding------------------------------------------------------------------
W.Ns1star=p(15);
W.a3star=p(16);
W.apf=D.kpf;
W.leo=D.leo;
W=pmac_winding(W,G);
c2=lte(W.dc*D.aso,G.wso); 

% Mass---------------------------------------------------------------------
mss = G.vsl*S.row;
mrs = G.vrb*R.row;
mpm = G.vpm*M.row*10;
mcd = 3*W.vcd*C.row;   
m   = mss+mrs+mpm+mcd;
c3  = lte(m,D.mlim);

% Circumscribing Volume----------------------------------------------------
tot_vol = (pi*G.rss^2)*(G.l+2*(W.leo+W.lew));

% Size---------------------------------------------------------------------
size = sqrt(m*tot_vol);                    % sqrt(kg m^3)


% check constraints satisified against those imposed-----------------------
CS = c1+c2+c3;
CI = 3;
if (CS<CI)
    f=eps*[1; 1]*(CS-NC)/NC;
    return
end

% controllaw---------------------------------------------------------------
Kt     = [p(17) p(18) p(19) p(20) p(21) p(22)];
Te_des = Kt.*D.Te; 
E  = pmac_electrical_parameters2(M,C,G,W);
wr = D.wrm*(G.P/2);
vll_max = D.vdc-2*D.vfs;
T = controllaw(Te_des,wr,G.P,E,vll_max);

% Lumped Parameters--------------------------------------------------------
vqsr  = E.Rs.*T.iqsr+(E.Ld.*T.idsr+E.lm).*wr;
vdsr  = E.Rs.*T.idsr-E.Lq.*T.iqsr.*wr;
vpkll = sqrt(3*(vqsr.^2+vdsr.^2));
Te = 0.75*G.P*(E.lm.*T.iqsr+(E.Ld-E.Lq).*T.iqsr.*T.idsr);

% check constraints satisified against those imposed-----------------------
sumc1 = 0;                                       
for i=1:6
    c4 = lte(T.idsr(i),0);
    c5 = lte(vpkll(i),vll_max);
    sumc1 = sumc1+c4+c5;
end
sumc1=sumc1/6;
CS = CS+sumc1;
CI = CI+2;
if (CS<CI)
    f=eps*[1; 1]*(CS-NC)/NC;
    return
end

% Current------------------------------------------------------------------
I.iqsr=T.iqsr;
I.idsr=T.idsr;
I=pmac_current(I);
I2345 = [I.Is(2) I.Is(3) I.Is(4) I.Is(5)];
Iopmax = max(I2345);
c6 = lte(Iopmax/W.ac,15*10^6);                % Current density constraint

% check constraints satisified against those imposed-----------------------
CS=CS+c6;
CI=CI+1;
if (CS<CI)
   f=eps*[1; 1]*(CS-NC)/NC;
   return
end

% Magnetic Field Analysis and Constraints----------------------------------
I0.Is=0;
I0.phii=0;
I0.iqsr=0;
I0.idsr=0;
sumc2 = 0;
for i =1:6
    F0 = pmac_field_analysis2(M,G,W,I0.Is,I0.phii,S,D.J,D.wrm(i));
    c7 = lte(max(abs(F0.Bt1c)),S.Blim);
    c8 = lte(max(abs(F0.Bb1c)),S.Blim);
    c9 = lte(F0.Brbtmx,R.Blim);
   c10 = lte(F0.Brbrmx,R.Blim);
   c11 = gte(F0.Hmn,M.Hci*D.km);
 sumc2 = sumc2+c7+c8+c9+c10+c11;
end
sumc2=sumc2/6;

% Check Constraints Satisified Against Those Imposed-----------------------
CS=CS+sumc2;
CI=CI+5;
if (CS<CI)
   f=eps*[1; 1]*(CS-NC)/NC;
   return
end

% Magnetic Field Analysis and Constraints----------------------------------
sumc3 = 0;
Tec = zeros(1,6);                               
Pr = zeros(1,6);                                
Ps = zeros(1,6);                                
Pc = zeros(1,6);                                
Pl = zeros(1,6);                                  
for i=1:6   
    F = pmac_field_analysis2(M,G,W,I.Is(i),I.phii(i),S,D.J,D.wrm(i));
 Pc(i)= F.Pc;
  c12 = lte(max(abs(F.Bt1c)),S.Blim);
  c13 = lte(max(abs(F.Bb1c)),S.Blim);
  c14 = lte(F.Brbtmx,R.Blim);
  c15 = lte(F.Brbrmx,R.Blim);
  c16 = gte(F.Hmn,M.Hci*D.km);

% Compute Corrected Torque-------------------------------------------------
if abs(D.wrm(i))>0
   Tec(i)= Te(i)-Pc(i)/D.wrm(i);
else
   Tec(i)= Te(i);
end
  c17 = gte(Tec(i),D.Te(i));
sumc3 = sumc3+c12+c13+c14+c15+c16+c17;
end
sumc3=sumc3/6;

% Check Constraints Satisified Against Those Imposed-----------------------
CS=CS+sumc3;
CI=CI+6;
if (CS<CI)
   f=eps*[1; 1]*(CS-NC)/NC;
   return
end

% Compute Total Loss------------------------------------------------------
sumc4 = 0;
for i=1:6
    Pr(i)  = 3.0*E.Rs*I.Is(i)^2;
    Ps(i)  = sixrt2opi*I.Is(i)*D.vfs;
    Pl(i)  = Pr(i)+Ps(i)+Pc(i);
    if i>=2 && i<=5                      % OP 1 and 6 transients therefore
        c18 = lte(Pl(i),D.Pllim);        % not included in loss constraint 
        sumc4 = sumc4+c18;
    end
end
sumc4=sumc4/4;

% Check Constraints Satisified Against Those Imposed-----------------------
CS=CS+sumc4;
CI=CI+1;
if (CS<CI)
   f=eps*[1; 1]*(CS-NC)/NC;
   return
end

% Compute Aggregate power--------------------------------------------------
Plagg = Pl(1)*0.05+Pl(2)*0.15+Pl(3)*0.05...         
        +Pl(4)*0.2+Pl(5)*0.4+Pl(6)*0.15;

Pin = zeros(1,6);
Pintotal = zeros(1,6);
Pout = zeros(1,6);
for i =1:6
    Pintotal(i) = (3/2)*(vqsr(i)*I.iqsr(i)+...
                    vdsr(i)*I.idsr(i))+Ps(i);
         Pin(i) = (3/2)*(vqsr(i)*I.iqsr(i)+...
                    vdsr(i)*I.idsr(i));
       Pout(i)  = Te(i)*D.wrm(i)-Pc(i);     
                
end
Pinagg =Pin(1)*0.05+Pin(2)*0.15+Pin(3)*0.05...        % Aggregate Input 
        +Pin(4)*0.2+Pin(5)*0.4+Pin(6)*0.15;
    
Pinaggtot =Pintotal(1)*0.05+Pintotal(2)*0.15...     
           +Pintotal(3)*0.05+Pintotal(4)*0.2+...
           Pintotal(5)*0.4+Pintotal(6)*0.15;
    
Poutagg = Pout(1)*0.05+Pout(2)*0.15+Pout(3)*0.05+...  
          Pout(4)*0.2+Pout(5)*0.4+Pout(6)*0.15;

efficiency = Poutagg/Pinaggtot;                     
eff_noinverter = Poutagg/Pinagg;

% Fitness------------------------------------------------------------------
f=[1/size; 1/Plagg];

% Reporting----------------------------------------------------------------
if (nargin==3)&&(fn>0)

   sd=3;
   disp([   ]);
   disp(['Design Data-------------------------------------------------']);
   disp(['    Outside Diameter: ' num2str(2*G.rss*100,sd) ' cm']);
   disp(['    Total Length: ' ...
              num2str((2*(W.leo+W.lew)+G.l)*100,sd) ' cm']);
   disp(['    Active Length: ' num2str(G.l*100,sd) ' cm']);
   disp(['    Number of Poles: ' num2str(G.P)]);
   disp(['    Number of Slots: ' num2str(G.Ss)]);
   disp(['    Stator Material Type: ' S.desc]);
   disp(['    Rotor Material Type: '  R.desc]);
   disp(['    Conductor Type: ' C.desc]);
   disp(['    Permanent Magnet Type: ' M.desc]);
   disp(['    Permanent Magnet Fraction: ' ...
              num2str(G.alphapm*100,sd) '%']);
   disp(['    Permanent Magnet Depth: ' ...
              num2str(G.dm*100,sd) ' cm']);
   disp(['    Shaft Radius: ' num2str(G.rrs*100,sd), ' cm']);
   disp(['    Inert Radius: ' num2str(G.rri*100,sd), ' cm']);
   disp(['    Rotor Iron Radius: ' num2str(G.rrb*100,sd) ' cm']);
   disp(['    Air Gap: ' num2str(G.g*1000,sd) ' mm']); 
   disp(['    Slot Depth: ' num2str((G.dtb+G.dttc)*100,sd) ' cm']);
   disp(['    Tooth Fraction: ' num2str(G.alphat*100,sd) ' %']);
   disp(['    Stator Backiron Depth: ' num2str(G.dsb*100,sd) ' cm']);
   disp(['    Rotor Backiron Depth: ' num2str(G.drb*100,sd) ' cm']);
   disp(['    Fundamental Conductor Density: ' ...
              num2str(W.Ns1,sd) ' cond/rad']);
   disp(['    3rd Harmonic Conductor Density: ' ....
              num2str(W.a3*100,sd) '%']);
   disp(['    Conductor Diameter: ' ...
              num2str(2*sqrt(W.ac/pi)*1000,sd), ' mm']);
   disp(['    Stator Iron Mass: ' num2str(mss,sd) ' kg']);
   disp(['    Rotor Iron Mass: ' num2str(mrs,sd) ' kg']);
   disp(['    Conductor Mass: ' num2str(mcd,sd) ' kg']);
   disp(['    Magnet Mass: ' num2str(mpm,sd) ' kg']);
   disp(['    Mass: '  num2str(m,sd) ' kg']);
   disp(['    A-Phase Winding Pattern (1st Pole): ']);
   disp(W.Nas(1:G.Ss/G.P));
   disp(['    Minimum Conductors Per Slot: ' num2str(min(W.Ns))]);
   disp(['    Maximum Conductors Per Slot: ' num2str(max(W.Ns))]);
   disp(['    Packing Factor: ' num2str(W.apf*100,sd) '%']);
   
   disp(['  ']); 
   disp(['Electrical Model--------------------------------------------']);
   disp(['    Number of Poles: ' num2str(G.P,sd)]);
   disp(['    Nominal Stator Resistance: ' ...
              num2str(E.Rs*1000,sd), ' mOhms']);
   disp(['    Q-Axis Inductance: ' num2str(E.Lq*1e3,sd) ' mH']);
   disp(['    D-Axis Inductance: ' num2str(E.Ld*1e3,sd) ' mH']);
   disp(['    Flux Linkage Due to PM: ' num2str(E.lm*1000,sd) ' mVs']);
      
   disp(['  ']);
   disp(['Operating Point Performance Data----------------------------']);
   disp(['    Speed: ' num2str(D.wrm*30/pi) ' RPM']); 
   disp(['    Frequency: ' num2str(D.wrm*(G.P/2)*(1/(2*pi)),sd) 'Hz']);
   disp(['    Q-Axis Voltage: ' num2str(vqsr,sd) ' V']);
   disp(['    D-Axis Voltage: ' num2str(vdsr,sd) ' V']);
   disp(['    Peak Line-to-Line Voltage: ' num2str(vpkll,sd) ' V']);
   disp(['    Q-Axis Current: ' num2str(I.iqsr,sd) ' A']);
   disp(['    D-Axis Current: ' num2str(I.idsr,sd) ' A']);
   disp(['    Peak Line Current: ' num2str(I.Is*sqrt(2),sd) ' A']);
   disp(['    Current Density: ' ...
              num2str(I.Is/(W.ac*1e6),sd) ' A rms/mm^2']);
   disp(['    Torque: ' num2str(Te,sd) ' Nm']);
   disp(['    Corrected torque: ' num2str(Tec,sd) ' Nm']);
   disp(['    Semiconductor Conduction Loss for each operating point: '...
              num2str(Ps,sd) ' W']);
   disp(['    Machine Resistive Losses for each operating point: '...
              num2str(Pr,sd) ' W']);
   disp(['    Machine Core Loss for each operating point: '...
              num2str(Pc,sd) ' W']);
   disp(['    Total Loss for each operating point: '...
              num2str(Pl,sd) ' W']);
   disp(['    Machine Aggregate Loss: ' ...
              num2str(Plagg,sd) 'W']);
   disp(['    Machine Aggregate Efficiency: ' ...
              num2str(eff_noinverter*100,sd) '%']);
   disp(['    Inverter Efficiency: ' ...
              num2str((Pinagg/Pinaggtot)*100,sd) '%']);
   disp(['    Machine/Inverter Efficiency: ' ...
              num2str(efficiency*100,sd) '%']);
   disp(['    Stator Tooth Flux Density / Limit: ' ...
              num2str(100*max(abs(F.Bt1c))/S.Blim,sd) '%']); 
   disp(['    Stator Backiron Flux Density / Limit: ' ...
              num2str(100*max(abs(F.Bb1c))/S.Blim,sd) '%']);
   disp(['    Rotor Peak Tangential Flux Density / Limit: ' ...
              num2str(100*F.Brbtmx/R.Blim,sd) '%']);
   disp(['    Rotor Peak Radial Flux Density / Limit: ' ...
              num2str(100*F.Brbrmx/R.Blim,sd) '%']);
   disp(['    Permanent Magnet Demagnetization / Limit: ' ...
              num2str(100*F.Hmn/(M.Hci*D.km),sd) '%']);
   
   % draw the machine
   pmac_plotting_v0_2(D,G,W,R,S,M,C,[fn fn+1 fn+2]);
   
   % draw field plots
   for i=1:6
       F = pmac_field_analysis2(M,G,W,I.Is(i),I.phii(i),S,D.J*100,D.wrm(i));
       figure(fn+3+i-1); 
       subplot(2,1,2);
       plot(F.qrc*180/pi,F.Bb1c,'g-',F.qrc*180/pi,F.Bt1c,'b-');
       legend('backiron','tooth');
       xlabel('\theta_r, Degrees');
       ylabel('Flux Density, T')
       set(gca,'XTick',[0 60 120 180 240 300 360]);
       a=axis;
       axis([0 360 a(3:4)]);
       grid on;
   end
       
  
end

% return design data in f
if (nargin==3)&&(fn==0)
   
   clear f;
   f.size =size;
   f.m=m;
   f.mss=mss;
   f.mrs=mrs;
   f.mpm=mpm;
   f.mcd=mcd;
   f.Pr1=Pr(1);
   f.Pr2=Pr(2);
   f.Pr3=Pr(3);
   f.Pr4=Pr(4);
   f.Pr5=Pr(5);
   f.Pr6=Pr(6);
   f.Ps1=Ps(1);
   f.Ps2=Ps(2);
   f.Ps3=Ps(3);
   f.Ps4=Ps(4);
   f.Ps5=Ps(5);
   f.Ps6=Ps(6); 
   f.Pc1=Pc(1);
   f.Pc2=Pc(2);
   f.Pc3=Pc(3);
   f.Pc4=Pc(4);
   f.Pc5=Pc(5);
   f.Pc6=Pc(6);
   f.Plagg=Plagg;
   f.P=G.P;
   f.Is1=I.Is(1);
   f.Is2=I.Is(2);
   f.Is3=I.Is(3);
   f.Is4=I.Is(4);
   f.Is5=I.Is(5);
   f.Is6=I.Is(6);
   f.J1=I.Is(1)/W.ac;
   f.J2=I.Is(2)/W.ac;
   f.J3=I.Is(3)/W.ac;
   f.J4=I.Is(4)/W.ac;
   f.J5=I.Is(5)/W.ac;
   f.J6=I.Is(6)/W.ac;
   f.ac=W.ac;
   f.st=p(1);
   f.rt=p(2);
   f.ct=p(3);
   f.mt=p(4);
   f.aslt=G.aslt;
   f.l=G.l;
   f.rss=G.rss;
   f.alphat=G.alphat;
   f.alphapm=G.alphapm;
   f.efficiency=efficiency;
   f.Pl1=Pl(1);
   f.Pl2=Pl(2);
   f.Pl3=Pl(3);
   f.Pl4=Pl(4);
   f.Pl5=Pl(5);
   f.Pl6=Pl(6);
   
end

end  % pmac_fitness

function [c]=lte(x,xmx)
%   FUNCTION:    lte(x,xmx)
%   DESCRIPTION: Function to test constraint. Returns 1 if the argument 
%                is less-than the allowed value, otherwise it returns 
%                a value between 0 and 1.
%   INPUTS:      x       -   variable to be measured
%                xmx     -   maximum allowed value for the variable
%   OUTPUTS:     f       -   fitness value
    if (x<=xmx)
        c=1;
    else
        c=1/(x-xmx+1);
    end
end

function [c]=gte(x,xmn)
%   FUNCTION:    [c]=gte(x,xmn)
%   DESCRIPTION: Function to test constraint. Returns 1 if the argument
%                is greater-than allowed value, otherwise it returns
%                a value between 0 and 1.
%   INPUTS:      x       -   variable to be measured
%                xmn     -   minimum allowed value for the variable
%   OUTPUTS:     f       -   fitness value

    if (x<xmn)
        c=1/(xmn-x+1);
    else
        c=1;
    end
end