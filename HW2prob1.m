% Modified By Avyay Sah
% Email: asah@purdue.edu

% Initialize
clear all; close all; clc

% Define some handy numbers
mm = 1e-3;
mu0 = 4*pi*10^-7;
cu_res =1.68*10^-8;      %resistivity of copper at room temp

% Specify the dimensions of the EI Core based on
% J.Cale, S.D.Sudhoff, and L.Tan, "Accurately Modeling EI Core Inductors 
% Using a High-Fidelity Magnetic Equivalent Circuit Approach", 
% IEEE Transactions on Magnetics, Vol.42, No.1, January 2006

wl = 25.0 * mm;
wdc = 60.3 * mm;
wd = 50.0 * mm;
ws = 30.0 * mm;
ds = 60.0 * mm;
g  = 5.00 * mm;
d  = 50.0 * mm;
c = 15.36 * mm;      %clearance

kpf = 0.6;                                        % packing factor
SA  = kpf * 40 * 40;                              % area of conductor(mm^2)                     
vol = 2*kpf*pi*40*mm*(((wd/2+c+40*mm)^2)-...      % volume of conductor
      ((wd/2+c)^2)) 
lcu = 2*pi*(wd/2+c+20*mm);
Rcu = cu_res*lcu/(SA*mm*mm);
% Curent density range of interest
j  = linspace(0,6,100);                         %current density in mm^2      
Npoints = length(j);

% Power Loss
PL = ((j*(10^6)).^2)*vol*cu_res;
% Determine the MMF source on the circuit
F1 = SA*j;
F2 = F1;

% Calculate the permeances on the air gap
Pfoc = mu0*d*log(1+pi*c/g)/pi;           
Pff = Pfoc;
Pgd = (mu0*wd*d)/g;
Pf = 2*(Pfoc+Pff);                          %total fringing
Pg = Pgd+Pf;                                %total air gap permeance 
Rgd = 1/Pgd;
Rf = 1/Pf;

%HIPERCO50 material parameters to calculate permeability as a function of B
HIPERCO50.mur = 43371.9609;
HIPERCO50.muB.a   = [0.43708   0.0003068  0.00026279  0.00024516];
HIPERCO50.muB.b   = [17.13367      2.139356      163.4348      1.476588];
HIPERCO50.muB.g   = [2.2836      1.3692      1.6772       3.494];
HIPERCO50.muB.d   = HIPERCO50.muB.a./HIPERCO50.muB.b;
HIPERCO50.muB.t   = exp(-HIPERCO50.muB.b.*HIPERCO50.muB.g);
HIPERCO50.muB.h   = HIPERCO50.muB.a.*HIPERCO50.muB.t;
HIPERCO50.muB.e   = HIPERCO50.muB.t./(HIPERCO50.muB.t+1);
HIPERCO50.muB.z   = 1./(HIPERCO50.muB.t+1);

% Initial flux linkage and convergence
B=NaN(Npoints,1);
c=NaN(Npoints,1);

% Initialize the MEC
MEC = mec_init(2,9,1);    
   
% Establish a list of non-linear materials
MEC = mec_nl_material(MEC,1,@muB,HIPERCO50);

% Branches
MEC= mec_nls_branch(MEC,1,wd*d,(ds/2+ws/2),1,0.0,0.0,[1 -2],0.0);  % MB 1
MEC= mec_lp_branch(MEC,2,Pg,[1 -2]);                               % MB 2
MEC= mec_nls_branch(MEC,3,wd*d,(ds/2+ws/2),1,0.0,0.0,[1 -2],0.0);  % MB 3
MEC= mec_nl_branch(MEC,4,ws*d,wd/2+wdc+wl/2,1,[1],0.0);            % MB 4
MEC= mec_nl_branch(MEC,5,wl*d,ws+2*ds+g,1,[1],0.0);                % MB 5
MEC= mec_nl_branch(MEC,6,ws*d,wd/2+wdc+wl/2,1,[1],0.0);            % MB 6
MEC= mec_nl_branch(MEC,7,ws*d,wd/2+wdc+wl/2,1,[2],0.0);            % MB 7
MEC= mec_nl_branch(MEC,8,wl*d,ws+2*ds+g,1,[2],0.0);                % MB 8
MEC= mec_nl_branch(MEC,9,ws*d,wd/2+wdc+wl/2,1,[2],0.0);            % MB 9 

% Loop over current values
for n = 1:Npoints

   % Assign MMF
   MEC.mb.Fs(1)=F1(n);
   MEC.mb.Fs(3)=F2(n);

   % Solve the MEC
   [Pm,Pb,PR,Fb,converge] = mec_solve(MEC);
   
   % Compute air gap flux density
   PHI_ag = -PR(2)*(Rf/(Rf+Rgd));
   B(n) = PHI_ag/(wd*d)
   
   % Copy convergence
   c(n)=converge;
    
end

% Generate the lambda-vs-current curve
if (min(c)<1)
   error('Solution Did Not Converge');
end

% Plot flux density over air gap versus current density
figure(1)
plot(j,B)
xlabel('Current Density, (A/mm^2)')
ylabel('Flux Density Over Air Gap, (T)')
grid on
axis([0 6 0 2.5]);

%Plot Power Loss versus flux desnity over air gap
figure(2)
plot(B,PL)
xlabel('Flux Density Over Air Gap, (T)')
ylabel('Power Loss, (Watts)')
grid on
axis([0 2.2 0 600]);