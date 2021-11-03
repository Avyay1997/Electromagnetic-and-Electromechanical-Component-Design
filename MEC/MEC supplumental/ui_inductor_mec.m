function [lambda,fe,c] = ui_inductor_mec(UI,W,i)
% ui_inductor_mec  Magnetic equivalent circuit (mec) analysis of a 
%                  ui-core inductor.  Based on class notes for EE695
%
% [lambda,f,c] = ui_inductor_mec(UI,W,i)
%
% Inputs:
% UI        = Stucture of ui core parameters
%  UI.g    = Air gap (m)
%  UI.ds   = Slot depth (m)
%  UI.ws   = Slot width (m)
%  UI.wi   = I-core width (m)
%  UI.we   = End post width (m)
%  UI.wb   = Base width (m)
%  UI.d    = Core depth (m)
%  UI.pmub = Pointer to function return permeability and its derivative
%  UI.MP   = Structure of material parameters used by EI.pmub
% W        = Structure of winding parameters
%  W.N     = Number of turns
%  W.dw    = Depth of winding (m)
%  W.ww    = Width of winding (m)
% i        = Vector of currents at which to perform analysis (A)
%
% Outputs:
% lambda   = Vector of flux linkages corressponding to i(Vs)
% fe       = Vector of force on i core (negative is attractive) (Nm)
% c        = Vector whose elements are 1 if the analysis converged and
%            zero otherwise
%
% Internal:
% fopi     = 4/pi
% topi     = 2/pi
% pio4     = pi/4
% mu0      = Permeability of free space
% s        = Size of input current
% N        = Length of input current vector
% c        = Convergence flag vector
% P34      = Permeance from node 3 to 4 {across air gap}(H)
% P56      = Permeance from node 5 to 6 {across air gap)(H)
% Phl      = Permeance of horizontal leakage (H)
% Phb      = Permeance of horizontal bypass path (H)
% Pvl      = Vertical leakage permeance (H)
% Pvb      = Vertical bypass permeance (H)
% wt       = Used to calculate face bypass permeance (m)
% Pfb      = Face bypass permeance (H)
% k1       = Used to calculate face leakage permeance (m)
% rt2      = sqrt(2)
% xmax     = Used to calulate face leakage permeance (m)
% Pfl      = Face leakage permeamce (H)
% Pf       = Total face permeance (H)
% AB2      = Area of branch 2 (m^2)
% AB5      = Area of branch 5 (m^2)
% AB6      = Area of branch 6 (m^2)
% AB7      = Area of branch 7 (m^2)
% AB8      = Area of branch 8 (m^2)
% AB10     = Area of branch 10 (m^2)
% AB11     = Area of branch 11 (m^2)
% lB2      = Length of branch 2 (m)
% lB5      = Length of branch 5 (m)
% lB6      = Length of branch 6 (m)
% lB7      = Length of branch 7 (m)
% lB8      = Length of branch 8 (m)
% lB10     = Length of branch 10 (m)
% lB11     = Length of branch 11 (m)
% phib2    = Flux into branch 2 (through base) (Wb)
% phib5    = Flux into branch 5 (through i) (Wb)
% phib7    = Flux into branch 7 (top of left leg) (Wb)
% phib8    = Flux into branch 8 (bottom of left leg) (Wb)
% phib10   = Flux into branch 10 (top of right leg) (Wb)
% phib11   = Flux into branch 11 (bottom of right leg) (Wb)
% g        = Air gap (limited if two small) (m)
% PPB12wrtg= Partial of permeance in branch 12 w.r.t. g (H/m)
% PPvlwrtg = Partial of Pvl w.r.t. g (H/m)
% Pvbwrtg  = Partial of Pvb w.r.t. g (H/m)
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
pio4=pi/4;
mu0=pi*4e-7;

% initalize flux linkage, convergence, and determine
% number of data points
s=size(i);
N=length(i);
lambda=zeros(s);
fe    =zeros(s);
c     =zeros(s);

% Compute air gap permeances
P34=mu0*(UI.we*UI.d/UI.g+...
         UI.d*log(1+pi*UI.wi/UI.g)/pi+...
         topi*UI.d*log(1+pio4*min(UI.ws,UI.ds)/UI.g)+...
         topi*UI.we*log(1+pi*UI.wi/UI.g));
P56=P34;

% Compute horizontal permeances
Phl=(mu0/3)*W.dw*(UI.d/UI.ws+topi*log(1+pi*UI.we/UI.ws));
Phb=mu0*(UI.ds-W.dw)*(UI.d/UI.ws + topi*log(1+pi*UI.we/UI.ws));

% Compute vertical permeances
Pvl=mu0*UI.d*W.ww/(12*(UI.ds+UI.g));
Pvb=0.5*mu0*(UI.ws-W.ww)*( UI.d/(UI.ds+UI.g)+ ...
                          topi*log(1+pi*min(UI.wb,UI.wi)/(UI.ds+UI.g)));

  
% Compute face permeance
wt=UI.we+(UI.ws-W.ww)/2.0;
Pfb=0.5*mu0*wt*(UI.d+2.0*UI.wb)/(W.dw+sqrt((wt+W.ww)*wt));
k1=abs(2.0*W.dw-W.ww);
rt2=sqrt(2.0);
xmax=rt2*min(0.5*W.ww,W.dw);
Pfl=mu0*(UI.d+2.0*UI.wb)/(16.0*W.ww^2*W.dw^2);
if (k1>0)
   Pfl=Pfl*(xmax^4+rt2*k1*xmax^3+0.25*(k1*xmax)^2- ...
            rt2*k1^3*xmax/8+(k1^4)*log(1+2*rt2*xmax/k1)/16);
else
   Pfl=Pfl*(xmax^4+rt2*k1*xmax^3+0.25*(k1*xmax)^2- ...
            rt2*(k1^3)*xmax/8);
end
Pf=Pfl+Pfb;

% Areas and length of branches
AB2=UI.wb*UI.d;
lB2=UI.ws+UI.we;
AB10=UI.we*UI.d;
lB10=0.5*(W.dw+UI.wb);
AB7=AB10;
lB7=lB10;
AB11=UI.we*UI.d;
lB11=UI.ds-0.5*W.dw;
AB8=AB11;
lB8=lB11;
AB5=UI.wi*UI.d;
lB5=lB2;

% Initial guesses on flux linkage
phib2=0;
phib5=0;
phib7=0;
phib8=0;
phib10=0;
phib11=0;

% Set up and solve for the flux linkage
for n=1:N
        
   % build and solve mec 
   [MEC] = mec_init(6,13);
   [MEC] = mec_lp_branch(  MEC, 1, Pf,[1]);
   [MEC] = mec_nls_branch( MEC, 2, AB2, lB2,UI.pmub,UI.MP,W.N*i(n),0,[2 -1],phib2);
   [MEC] = mec_lp_branch(  MEC, 3, Phl,[3 -2]);
   [MEC] = mec_lp_branch(  MEC, 4, Phb,[4 -3]);
   [MEC] = mec_nl_branch(  MEC, 5, AB5, lB5,UI.pmub,UI.MP,[4],phib5);
   [MEC] = mec_lp_branch(  MEC, 6, Pvb,[5]);
   [MEC] = mec_nl_branch(  MEC, 7, AB7, lB7,UI.pmub,UI.MP,[2 -5],phib7);
   [MEC] = mec_nl_branch(  MEC, 8, AB8, lB8,UI.pmub,UI.MP,[3 -5],phib8);
   [MEC] = mec_lp_branch(  MEC, 9, P56, [4 -5]);
   [MEC] = mec_nl_branch(  MEC,10,AB10,lB10,UI.pmub,UI.MP,[6 -2],phib10);
   [MEC] = mec_nl_branch(  MEC,11,AB11,lB11,UI.pmub,UI.MP,[6 -3],phib11);
   [MEC] = mec_lp_branch(  MEC,12, P34,[6 -4]);
   [MEC] = mec_lp_branch(  MEC,13, Pvb,[6]);
  
   [Pm,Pb,PR,Fb,converge] = mec_solve(MEC);
   
   % solution for present current is used as initial condition
   % for next solution
   phib2=PR(2);
   phib5=PR(5);
   phib7=PR(7);
   phib8=PR(8);
   phib10=PR(10);
   phib11=PR(11);
      
   % compute flux linkage
   lambda(n)=W.N*(-PR(2)+W.N*i(n)*Pvl);
   
   % copy convergence
   c(n)=converge;
        
   % compute force
   g=UI.g;
   if (g<1e-9)
      g=1.e-9;
   end
   PPB12wrtg=(-mu0/g^2)* ...
             (UI.we*UI.d + ...
              UI.d*g*UI.wi/(g+pi*UI.wi)+ ...
              2*UI.d*g*min(UI.ws,UI.ds)/(4*g+pi*min(UI.ws,UI.ds))+ ...
              2.0*UI.we*UI.wi*g/(g+pi*UI.wi));
   PPvlwrtg=-mu0*UI.d*W.ww/(12*(UI.ds+g)^2); 
   PPvbwrtg=-(0.5*mu0*(UI.ws-W.ww)/(UI.ds+g)^2)* ...
             (UI.d+2.0*(UI.ds+g)*min(UI.wb,UI.wi)/(UI.ds+g+pi*min(UI.wb,UI.wi)));
   fe(n)=PPB12wrtg*Fb(12)^2+PPvbwrtg*Fb(6)^2+PPvlwrtg*(W.N*i(n))^2;

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
