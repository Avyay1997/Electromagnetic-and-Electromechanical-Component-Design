function [Bg,c] = EE_mec(EE,W,J)
% EE_mec   Magnetic equivalent circuit (mec) analysis of a 
%                  ei-core inductor.  Based on paper:
%
%                  J.L. Cale, S.D. Sudhoff, “Accurately Modeling EI Core 
%                  Inductors Using a High Fidelity Magnetic Equivalent
%                  Circuit Approach,” IEEE Transactions on Magnetics, 
%                  Vol. 42, Issue 1, Jan 2006, pp. 40-46.
%
%                  and class notes for EE695
%
% [Bg,c] = EE_mec(EE,W,J)
%

% Modified By Avyay Sah
% Email: asah@purdue.edu

% Define some handy numbers
mu0 = 4*pi*10^-7;

% Determine the MMF source on the circuit
F1 = W.kpf*W.dw*W.ww*J;
F2 = F1;

% Calculate the permeances on the air gap
Pfoc = mu0*EE.lc*log(1+pi*(EE.ds/2)/EE.g)/pi;           
Pff = Pfoc;
Pgd = (mu0*EE.wc*EE.lc)/EE.g;
Pf = 2*(Pfoc+Pff);                          %total fringing
Pg = Pgd+Pf;                                %total air gap permeance 
Rgd = 1/Pgd;
Rf = 1/Pf;

% Initial flux linkage and convergence
Npoints=length(J);
Bg=NaN(Npoints,1);
c=NaN(Npoints,1);

% Initialize the MEC
MEC = mec_init(2,9,1);    
   
% Establish a list of non-linear materials
MEC = mec_nl_material(MEC,1,@muB,EE.MP);

% Branches
MEC= mec_nls_branch(MEC,1,EE.wc*EE.lc,(EE.ds/2+EE.wb/2),1,0.0,0.0,[1 -2],0.0);     % MB 1
MEC= mec_lp_branch(MEC,2,Pg,[1 -2]);                                               % MB 2
MEC= mec_nls_branch(MEC,3,EE.wc*EE.lc,(EE.ds/2+EE.wb/2),1,0.0,0.0,[1 -2],0.0);     % MB 3
MEC= mec_nl_branch(MEC,4,EE.wb*EE.lc,EE.wc/2+EE.ws+EE.we/2,1,[1],0.0);             % MB 4
MEC= mec_nl_branch(MEC,5,EE.we*EE.lc,EE.wb+2*EE.ds+EE.g,1,[1],0.0);                % MB 5
MEC= mec_nl_branch(MEC,6,EE.wb*EE.lc,EE.wc/2+EE.ws+EE.we/2,1,[1],0.0);             % MB 6
MEC= mec_nl_branch(MEC,7,EE.wb*EE.lc,EE.wc/2+EE.ws+EE.we/2,1,[2],0.0);             % MB 7
MEC= mec_nl_branch(MEC,8,EE.we*EE.lc,EE.wb+2*EE.ds+EE.g,1,[2],0.0);                % MB 8
MEC= mec_nl_branch(MEC,9,EE.wb*EE.lc,EE.wc/2+EE.ws+EE.we/2,1,[2],0.0);             % MB 9 

% Loop over current values
for n = 1:Npoints

   % Assign MMF
   MEC.mb.Fs(1)=F1(n);
   MEC.mb.Fs(3)=F2(n);

   % Solve the MEC
   [Pm,Pb,PR,Fb,converge] = mec_solve(MEC);
   
   % Compute air gap flux density
   PHI_ag = -PR(2)*(Rf/(Rf+Rgd));
   Bg(n) = PHI_ag/(EE.wc*EE.lc);
   
   % Copy convergence
   c(n)=converge;
    
end
    
end
