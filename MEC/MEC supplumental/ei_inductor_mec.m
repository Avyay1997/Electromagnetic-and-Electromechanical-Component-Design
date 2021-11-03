function [lambda,c] = ei_inductor_mec(EI,W,i)
% ei_inductor_mec  Magnetic equivalent circuit (mec) analysis of a 
%                  ei-core inductor.  Based on paper:
%
%                  J.L. Cale, S.D. Sudhoff, “Accurately Modeling EI Core 
%                  Inductors Using a High Fidelity Magnetic Equivalent
%                  Circuit Approach,” IEEE Transactions on Magnetics, 
%                  Vol. 42, Issue 1, Jan 2006, pp. 40-46.
%
%                  and class notes for EE695
%
% [lambda,c] = ei_inductor_mec(EI,W,i)
%
% Inputs:
% EI       = Stucture of ei core inductor parameters
%  EI.g    = Air gap (m)
%  EI.ds   = Slot depth (m)
%  EI.ws   = Slot width (m)
%  EI.wi   = I-core width (m)
%  EI.wc   = Center post width (m)
%  EI.we   = End post width (m)
%  EI.wu   = Underneeth section width (m)
%  EI.d    = Core depth (m)
%  EI.pmub = Pointer to function return permeability and its derivative
%  EI.MP   = Structure of material parameters used by EI.pmub
% W        = Structure of winding parameters
%  W.N     = Number of turns
%  W.dw    = Depth of winding (m)
%  W.ww    = Width of winding (m)
% i        = Vector of currents at which to perform analysis (A)
%
% Outputs:
% lambda   = Vector of flux linkages corressponding to i(Vs)
% c        = Vector whose elements are 1 if the analysis converged and
%            zero otherwise
%
% Internal:
% fopi     = 4/pi
% topi     = 2/pi
% mu0      = Permeability of free space
% s        = Size of input current
% N        = Length of input current vector
% c        = Convergence flag vector
% P23      = Permeance from node 2 to 3 {center leg to I-core)(H)
% P45      = Permeance from node 4 to 5 {end leg to I-core)(H)
% TP45     = 2*P45
% Phlp1    = Permeance of horizontal leakage path 1 (H)
% Phlp2    = Permeance of horizontal leakage path 2 (H)
% Phl      = Permeance of horizontal leakage (H)
% TPhl     = Twice permeance of horizontal leakage
% Phb      = Permeance of horizontal bypass path (H)
% TPhb     = Twice Permeance of horizontal bypass path (H)
% Pvl      = Vertical leakage permeance (H)
% TPvl     = Twice vertical leakage permeance (H)
% Pvb      = Vertical bypass permeance (H)
% TPvb     = Twice vertical bypass permeance (H)
% K        = Used to calculate face permeance
% xmax     = Used to calculate face permeance
% Pfl      = Face leakage permeance (internal) (H)
% wt       = used to calcuate find external face permeance (m)
% Rfb      = face bypass reluctance (1/H)
% Pf       = Total face permeance (H)
% TPf      = Twice face permeance (H)
% A01      = Area of branch from node 0 to 1 (bottom of center leg) (m^2)
% l01      = Lenght of branch from node 0 to 1 (bottom of center leg) (m)
% A12      = Area of branch from node 1 to 2 (top of center leg) (m^2)
% l12      = Length of branch from node 1 to 2 (top of center leg) (m)
% TA34     = Twice area of branch from node 3 to 4 (top of I) (m^2) 
% l34      = Lenght of branch from node 3 to 4 (m)
% TA56     = Twice area of branch from node 5 to 6 (top of outer leg) (m^2)
% l56      = Length of branch from node 5 to 6 (top of outer leg) (m)
% TA67     = Twice area from node 6 to 7 (bottom of outer leg) (m)
% l67      = Length of branch from node 6 to 7 (bottom of oter leg) (m)
% phib2    = Flux into branch 2 (bottom of center leg) (Wb)
% phib3    = Flux into branch 3 (top of center leg) (Wb)
% phib5    = Flux into branch 5 (I core in both directions) (Wb)
% phib8    = Flux into branch 8 (bottom of E in both direction) (Wb)
% phib10   = Flux into branch 10 (top of both outer legs) (Wb)
% pbhb11   = Flux into branch 11 (bottom of both outer legs) (Wb)
% n        = Number of current points
% MEC      = MEC structure
% Pm       = Vector of mesh fluxes (Wb)
% Pb       = Vector of branch fluxes (Wb)
% PR       = Vector of flux through reluctance branch element (Wb)
% Fb       = Vector of branch MMFs (A)
%
% Written by:
% S.D. Sudhoff
% School of Electrical and Computer Engineering
% 1285 Electrical Engineering Building
% West Lafayette, IN 47907-1285
% Phone: 765-494-3246
% Fax: 765-494-0676
% E-mail: sudhoff@ecn.purdue.edu

% constants
fopi=4.0/pi;
topi=2.0/pi;
mu0=pi*4e-7;

% initalize flux linkage, convergence, and determine
% number of data points
s=size(i);
N=length(i);
lambda=zeros(s);
c     =zeros(s);

% Compute air gap permeances
P23=mu0*(EI.d*EI.wc/EI.g + ...
         EI.d*fopi*log(1+0.25*pi*min(EI.ws,EI.ds)/EI.g) + ...
         EI.wc*topi*log(1+pi*EI.wi/EI.g));
P45=mu0*(EI.we*EI.d/EI.g+ ...
        (EI.d+2.0*EI.we)*log(1+pi*EI.wi/EI.g)/pi+ ...
         topi*EI.d*log(1+0.25*pi*min(EI.ws,EI.ds)/EI.g));
TP45=2*P45; 

% Compute horizontal bypass permeances
Phlp1=mu0*W.dw*EI.d/(3.0*EI.ws);
Phlp2=2.0*W.dw*mu0*log(1+pi*min(EI.we,0.5*EI.wc)/EI.ws)/(3.0*pi);
Phl=Phlp1+Phlp2;
TPhl=2.0*Phl;
Phb=mu0*EI.d*(EI.ds-W.dw)/EI.ws+ ...
    topi*mu0*(EI.ds-W.dw)*log(1+pi*min(EI.we,0.5*EI.wc)/EI.ws);
TPhb=2.0*Phb;

% Compute vertical leakage permeances
Pvl=mu0*W.ww*EI.d/(12*(EI.ds+EI.g));
TPvl=2.0*Pvl;
Pvb=mu0*EI.d*(EI.ws-W.ww)/(EI.ds+EI.g);
TPvb=2.0*Pvb;

% Compute face bybass reluctance
K=abs(W.dw-2.0*W.ww);
xmax=min(W.ww,0.5*W.dw)*sqrt(2.0);
Pfl=mu0*EI.wc*(xmax^4+sqrt(2)*K*xmax^3+0.25*xmax^2*K^2 - ...
               sqrt(2)*K^3*xmax/8+...
               K^4*log(1+2.0*sqrt(2.0)*xmax/K)/16.0)/ ...
               (16.0*W.ww^2*W.dw^2);
wt=EI.ds-W.dw;
if wt>0
   Rfb=(wt+EI.wu)*(W.ww+sqrt((2*EI.d+wt+EI.wu)*wt*EI.wu/(wt+EI.wu)))/ ...
       (mu0*wt*EI.wu*EI.wc);
   Pf=Pfl+1/Rfb;
else
   Pf=Pfl;
end
TPf=2.0*Pf;

% Areas and length of branches
A01=EI.wc*EI.d; 
l01=0.5*(EI.ds+EI.wu);

A12=EI.wc*EI.d;
l12=0.5*EI.ds;

TA34=2.0*EI.wi*EI.d;
l34=EI.ws+0.5*(EI.wc+EI.we);

TA07=2.0*EI.wu*EI.d;
l07=EI.ws+0.5*(EI.wc+EI.we);

TA56=2.0*EI.we*EI.d;
l56=0.5*EI.ds;

TA67=2.0*EI.we*EI.d;
l67=0.5*(EI.ds+EI.wu);

% Initial guesses on flux linkage
phib2=0;
phib3=0;
phib5=0;
phib8=0;
phib10=0;
phib11=0;

% Set up and solve for the flux linkage
for n=1:N
        
   % build and solve mec 
   [MEC] = mec_init(5,12);
   [MEC] = mec_lp_branch( MEC, 1, TPf,[1]);
   [MEC] = mec_nls_branch(MEC, 2, A01,l01,EI.pmub,EI.MP,W.N*i(n),0,[1 -4],phib2);
   [MEC] = mec_nl_branch( MEC, 3, A12,l12,EI.pmub,EI.MP,[1 -3],phib3);
   [MEC] = mec_lp_branch( MEC, 4, P23,[2]);
   [MEC] = mec_nl_branch( MEC, 5,TA34,l34,EI.pmub,EI.MP,[2],phib5);
   [MEC] = mec_lp_branch( MEC, 6,TPhb,[2 -3]);
   [MEC] = mec_lp_branch( MEC, 7,TPhl,[3 -4]);
   [MEC] = mec_nl_branch( MEC, 8,TA07,l07,EI.pmub,EI.MP,[4],phib8);
   [MEC] = mec_lp_branch( MEC, 9,TP45,[2 -5]);
   [MEC] = mec_nl_branch( MEC,10,TA56,l56,EI.pmub,EI.MP,[3 -5],phib10);
   [MEC] = mec_nl_branch( MEC,11,TA67,l67,EI.pmub,EI.MP,[4 -5],phib11);
   [MEC] = mec_lp_branch( MEC,12,TPvb,[5]);
   
   [Pm,Pb,PR,Fb,converge] = mec_solve(MEC);
   
   % solution for present current is used as initial condition
   % for next solution
   phib2=PR(2);
   phib3=PR(3);
   phib5=PR(5);
   phib8=PR(8);
   phib10=PR(10);
   phib11=PR(11);
   
   % compute flux linkage
   lambda(n)=W.N*(-PR(2)+W.N*i(n)*TPvl);
   
   % copy convergence
   c(n)=converge;
        
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
%  License along with MEC Toolbox.  If not, see
%  <http://www.gnu.org/licenses/>.
