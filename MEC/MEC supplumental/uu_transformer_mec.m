function [lambdap,lambdas,Bulp,Buls,Bub,c] = ...
          uu_transformer_mec(U,PW,SW,ip,is)
% uu_transformer_mec  Magnetic equivalent circuit (mec) analysis of a 
%                     uu-core transformer.  Note this code can also be
%                     used to analyze a UI transformer if the depth
%                     of the slot is taken to be 1/2 of the depth of the
%                     actual core. It is assumed the primary and secondary
%                     winding are each wrapped around one leg.
%
% [lambdap,lambdas,Bulp,Buls,Bub,c] = uu_transformer_mec(U,PW,SW,ip,is)
%
% Inputs:
% U        = stucture of uu core inductor parameters
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
%  U.MP    = structure of core magnetic material parameters
% PW       = structure of primary winding parameters
%  PW.lcr  = length of core (m)
%  PW.cl   = clearance on each side of core length (m)
%  PW.wcr  = width of core (m)
%  PW.cw   = clearance on each side of core width (m)
%  PW.dcr  = depth of core portion surrounded by core (m)
%  PW.cd   = clearance in the depth direction (m)
%  PW.wcl  = width of coil cross section (m)
%  PW.dcl  = depth of coil cross section (m)
%  PW.Nt   = number of turns
%  PW.acnd = area of conductor (m^2)
%  PW.lcnd = length of conductor (m)
%  PW.vcnd = volume of conductor (m^3)
%  PW.pf   = packing factor
%  PW.wwn  = width of window {defined as rectangle which is at least
%            C.cwn from the coil} (m)
%  PW.lwn  = length of window (m)
%  PW.cwn  = clearance to window (m)
%  PW.lex  = mean excess length of bundle (beyond perimeter of core)
%  PW.vol  = volume of coil (including dead space)
%  PW.gvld = 1 if valid geometry 0<=C.gvld<1 if not
% SW       = structure of secondary winding parameters
%  SW.lcr  = length of core (m)
%  SW.cl   = clearance on each side of core length (m)
%  SW.wcr  = width of core (m)
%  SW.cw   = clearance on each side of core width (m)
%  SW.dcr  = depth of core portion surrounded by core (m)
%  SW.cd   = clearance in the depth direction (m)
%  SW.wcl  = width of coil cross section (m)
%  SW.dcl  = depth of coil cross section (m)
%  SW.Nt   = number of turns
%  SW.acnd = area of conductor (m^2)
%  SW.lcnd = length of conductor (m)
%  SW.vcnd = volume of conductor (m^3)
%  SW.pf   = packing factor
%  SW.wwn  = width of window {defined as rectangle which is at least
%            C.cwn from the coil} (m)
%  SW.lwn  = length of window (m)
%  SW.cwn  = clearance to window (m)
%  SW.lex  = mean excess length of bundle (beyond perimeter of core)
%  SW.vol  = volume of coil (including dead space)
%  SW.gvld = 1 if valid geometry 0<=C.gvld<1 if not
% ip       = vector of primary currents at which to perform analysis (A)
% is       = vector of secondary currents at which to perform analysis (A)
%            Note each element of the secondary current vector corresponds
%            to an element of the primary current vector
%
% Outputs:
% lambdap  = vector of primary flux linkages corresponding to ip & is(Vs)
% lambdas  = vector of secondary flux linkges corresponding to ip & is (Vs)
% Bulp     = vector of flux density in u-core leg on primary side (T)
% Buls     = vector of flux density in u-core leg on secondary side (T)
% Bub      = vector of flux density in u-core base (T)
% c        = vector whose elements are 1 if the analysis converged and
%            zero otherwise
%
% Internal:
% Ppil1    = internal leakage permeance associated with outside 
%            portion of primary winding {along length of U}(H)
% Ppel1    = exteral leakage permeance associated with outside
%            portion of primary winding {along lengh of U} (H)
% Psil1    = internal leakage permeance associated with outside 
%            portion of secondary winding {along length of U}(H
% Psel1    = external leakage permeance associated with outside
%            portion of secondary winding {along lengh of U} (H)
% Ppil2    = internal leakage permeance associated with outside 
%            portion of primary winding {along faces of U}(H)
% Ppel2    = exteral leakage permeance associated with outside
%            portion of primary winding {along faces of U} (H)
% Psil2    = internal leakage permeance associated with outside 
%            portion of secondary winding {along faces of U}(H
% Psel2    = external leakage permeance associated with outside
%            portion of secondary winding {along faces of U} (H)
% Ppil3    = intenal leakage permeance associate with inside
%            portion of primary winding {inner face of U} (H)
% Psil3    = intenal leakage permeance associate with inside
%            portion of primary winding {inner face of U} (H)
% Ppil     = total primary winding internal leakage permeance (H)
% Ppel     = total primary winding external leakage permeance (H)
% Psil     = total secondary winding internal leakage permeance (H)
% Psel     = total secondary winding external leakage permeance (H)
% phiest2  = initial estimated branch 2 flux (Wb)
% phiest3  = initial estimated branch 3 flux (Wb)
% phiest4  = initial estimated branch 4 flux (Wb)
% Auuleg   = cross sectional area of a leg of UU core (m^2)
% luuleg   = lenght of a UU leg {include both cores} (m)
% Auubase  = cross sectional areal of U core base (m^2)
% luubase  = lenght of sume of UU base segements (m)
% htb_pw   = height from top of pri. winding bundle to top of U (m);
% htb_sw    =height from top of sec. winding bundle to top of U (m);
%
% Written by:
% S.D. Sudhoff
% School of Electrical and Computer Engineering
% 1285 Electrical Engineering Building
% West Lafayette, IN 47907-1285
% Phone: 765-494-3246
% Fax: 765-494-0676
% E-mail: sudhoff@ecn.purdue.edu

% initalize flux linkage, convergence, and determine
% number of data points
s=size(ip);
N=length(ip);
lambdap=zeros(s);
lambdas=zeros(s);
Bulp   =zeros(s); 
Buls   =zeros(s);
Bub    =zeros(s);
c      =zeros(s);

% Compute the internal and external leakage permeances for outside of
% the U portions of the winding
htb_pw=U.du+U.ds-0.5*PW.dcl;
htb_sw=U.du+U.ds-0.5*SW.dcl;
[Ppil1,Ppel1] = P_rct_bndle_fm(htb_pw,htb_pw, ...
                               PW.dcl,PW.wcl, ...
                               PW.cw,U.l);
[Psil1,Psel1] = P_rct_bndle_fm(htb_sw,htb_sw, ...
                               SW.dcl,SW.wcl, ...
                               SW.cw,U.l);
                           
% Compute the internal and external leakage permeances for faces of
% the U portions of the winding                           
[Ppil2,Ppel2] = P_rct_bndle_fm(htb_pw,htb_pw, ...
                               PW.dcl,PW.wcl, ...
                               PW.cl,2*U.wu);
[Psil2,Psel2] = P_rct_bndle_fm(htb_sw,htb_sw, ...
                               SW.dcl,SW.wcl, ...
                               SW.cl,2*U.wu);

% Compute the internal leakage permeance for inside the U and corner
% portions of winding
Ppil3 = Pint_rct_bndle_air(PW.dcl,PW.wcl,PW.lex);
Psil3 = Pint_rct_bndle_air(SW.dcl,SW.wcl,SW.lex);

% Compute total internal permeances
Ppil=Ppil1+Ppil2+Ppil3;
Psil=Psil1+Psil2+Psil3;
Ppel=Ppel1+Ppel2;
Psel=Psel1+Psel2;

% Areas and length of branches
Auuleg=U.l*U.wu;
luuleg=2.0*U.ds+U.du;
Auubase=U.l*U.du;
luubase=2.0*(U.ws+U.wu);

% Initial guesses on flux linkage
phiest2=0;
phiest3=0;
phiest4=0;

% Set up and solve for the flux linkage
for n=1:N
        
   % build and solve mec 
   [MEC] = mec_init(3,5);
   [MEC] = mec_lp_branch(MEC,1,Ppel,[-1]);
   [MEC] = mec_nls_branch(MEC,2,Auuleg, luuleg,U.pmub, ...
                          U.MP,PW.Nt*ip(n),0,[1 -2],phiest2);
   [MEC] = mec_nl_branch(MEC,3,Auubase,luubase,U.pmub,U.MP,[2],phiest3);
   [MEC] = mec_nls_branch(MEC,4,Auuleg, luuleg,U.pmub, ...
                          U.MP,SW.Nt*is(n),0,[2 -3],phiest4);
   [MEC] = mec_lp_branch(MEC,5,Psel,[3]);
   [Pm,Pb,PR,Fb,converge] = mec_solve(MEC);
   
   % solution for present current is used as initial condition
   % for next solution
   phiest2=PR(2);
   phiest3=PR(3);
   phiest4=PR(4);
   
   % compute flux linkage
   lambdap(n)=PW.Nt*(-PR(2)+PW.Nt*ip(n)*Ppil);
   lambdas(n)=SW.Nt*(-PR(4)+SW.Nt*is(n)*Psil);
   
   % compute flux densities
   Bub(n)=  PR(3)/(U.du*U.l);
   Bulp(n)=-PR(2)/(U.wu*U.l);
   Buls(n)= PR(4)/(U.wu*U.l);
      
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
