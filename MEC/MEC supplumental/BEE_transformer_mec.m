function [lambdap,lambdas,Bcp,Bop,Bb,c] = ...
         BEE_transformer_mec(BE,PBW,SBW,g,ip,is)
% BEE_transformer_mec  Magnetic equivalent circuit (mec) analysis of a 
%                      bobbin EE-core transformer.  Note this code can also be
%                      used to analyze a UI transformer if the depth
%                      of the slot is taken to be 1/2 of the depth of the
%                      actual core. It is assumed the primary and secondary
%                      winding are each wrapped around one leg.
%
% [lambdap,lambdas,Bcp,Bop,Bb,c] = BEE_transformer_mec(BE,PBW,SBW,ip,is)
%
% Inputs:
% BE        = structure of bobin E-core parameters
%  BE.rh    = radius of central hole (m)
%  BE.wip   = width of center post (m)
%  BE.ws    = width of slot (m)
%  BE.wop   = width of outer post (m)
%  BE.hp    = height of post (m)
%  BE.hb    = height of base (m)
%  BE.rcp   = outside radius of center post (m)
%  BE.rop   = inside radius of outers posts (m)
%  BE.re    = smallest radius to outside edge (m)
%  BE.qic   = angle from E centerline to inner corner of outer post(rad);
%  BE.qoc   = angle from E centerline to outer corner of outer post(rad);
%  BE.mcs_cp= magnetic cross section of center post (m^2);
%  BE.mcs_op= magnetic cross section of one edge post (m^2);
%  BE.mcs_b = magentic cross section of base (m^2);
%  BE.v     = volume of core (m^3)
%  BE.MP    = structure of material parameters (if MP supplied)
%  BE.m     = mass of core (kb) (if MP supplied)
%  BE.pmub  = pointer to permeability function (if supplied)
% PBW       = structure of primary bobin winding parameters
%  PBW.ri   = inner radius of bobbin (m)
%  PBW.ro   = outer radius of bobbin (m)
%  PBW.w    = width of bobbin cross section (m)
%  PBW.l    = mean length of bobbin (circumference of center) (m)
%  PBW.h    = height of bobbin (m)
%  PBW.N    = number of turns
%  PBW.ac   = conductor area (m^2)
%  PBW.aw   = wire area (m^2)
%  PBW.bcs  = bobbin cross section (m^2)
%  PBW.wpf  = wire packing factor 
%  PBW.cpf  = conductor packing factor
%  PBW.vc   = volume of conductor (m^3)
%  PBW.lc   = length of conductor (m^3)  
%  PBW.mc   = mass of conductor (m^3) {if CP supplied }
%  PBW.rdc  = dc resistance of conductor (Ohms) {if CP supplied}
% SBW       = structure of secondary bobin winding parameters
%  SBW.ri   = inner radius of bobbin (m)
%  SBW.ro   = outer radius of bobbin (m)
%  SBW.w    = width of bobbin cross section (m)
%  SBW.l    = mean length of bobbin (circumference of center) (m)
%  SBW.h    = height of bobbin (m)
%  SBW.N    = number of turns
%  SBW.ac   = conductor area (m^2)
%  SBW.aw   = wire area (m^2)
%  SBW.bcs  = bobbin cross section (m^2)
%  SBW.wpf  = wire packing factor 
%  SBW.cpf  = conductor packing factor
%  SBW.vc   = volume of conductor (m^3)
%  SBW.lc   = length of conductor (m^3)  
%  SBW.mc   = mass of conductor (m^3) {if CP supplied }
%  SBW.rdc  = dc resistance of conductor (Ohms) {if CP supplied}
% g         = airgap
% ip        = vector of primary currents at which to perform analysis (A)
% is        = vector of secondary currents at which to perform analysis (A)
%             Note each element of the secondary current vector corresponds
%             to an element of the primary current vector
%
% Outputs:
% lambdap  = vector of primary flux linkages corresponding to ip & is(Vs)
% lambdas  = vector of secondary flux linkges corresponding to ip & is (Vs)
% Bcp      = vector of flux density in center post (T)
% Bop      = vector of flux density in outer post (T)
% Bb       = vector of flux density in u-core base (T)
% c        = vector whose elements are 1 if the analysis converged and
%            zero otherwise
%
% Internal:
% Ppil     = total primary winding internal leakage permeance (H)
% Psil     = total secondary winding internal leakage permeance (H)
% Acpost
% lcpost
% Aopost
% lopost
% Abase
% lbase
% F
% phiest
% Pm,Pb,PR,Fb,converge
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
Bcp    =zeros(s); 
Bop    =zeros(s);
Bb     =zeros(s);
c      =zeros(s);

% Compute the internal leakage permeance 
Ppil = Pint_rct_bndle_air(PBW.h,PBW.w,PBW.l);
Psil = Pint_rct_bndle_air(SBW.h,SBW.w,SBW.l);

% Areas and length of branches
Acpost=BE.mcs_cp;
lcpost=2*BE.hp+BE.hb;
Aopost=BE.mcs_op*2;               
lopost=lcpost;
Abase=BE.mcs_b*2;                
lbase=2*(0.5*(BE.re+BE.rop)-0.5*BE.rcp);

% Air gap permeances
% mu0=1e-7*pi;
mu0=4e-7*pi;
rmax1=BE.hp/2;
rmax2=min(rmax1,BE.rh/2);
Pgcpost=mu0*(Acpost/g+ ...                       % direct path
             2.0*BE.rcp*log(1+pi*rmax1/g)+ ...   % fring on outside of post
             2.0*BE.rh*log(1+pi*rmax2/g));       % fringe on inside of post
Pgopost=mu0*(Aopost/g+...                        % direct path (both sides)
             2.0*(4*BE.rcp+2*BE.wop)*log(1+pi*rmax1/g)/pi); % fringe

% Initial guesses on flux linkage
phiest=0;

% Set up and solve for the flux linkage
for n=1:N
    
   % build and solve mec 
   [MEC] = mec_init(1,5);
   F=PBW.N*ip(n)-SBW.N*is(n);
   [MEC] = mec_nls_branch(MEC,1,Acpost, lcpost,BE.pmub, ...
                          BE.MP,F,0,1,phiest);
   [MEC] = mec_lp_branch(MEC,2,Pgcpost,1);
   [MEC] = mec_nl_branch(MEC,3,Aopost,lopost,BE.pmub,BE.MP,1,phiest);
   [MEC] = mec_lp_branch(MEC,4,Pgopost,1);
   [MEC] = mec_nl_branch(MEC,5,Abase , lbase,BE.pmub,BE.MP,1,phiest);
   [Pm,Pb,PR,Fb,converge] = mec_solve(MEC);
   
   % solution for present current is used as initial condition
   % for next solution
   phiest=PR(1);
   
   % compute flux linkage
   lambdap(n)=PBW.N*(-PR(1)+PBW.N*ip(n)*Ppil);
   lambdas(n)=SBW.N*( PR(1)+SBW.N*is(n)*Psil);
   
   % compute flux densities
   Bcp(n)= -PR(1)/Acpost;
   Bop(n)= -PR(1)/Aopost;
   Bb(n) = -PR(1)/Abase;
      
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
%  License along with MEC Toolbox.  If not, see <http://www.gnu.org/licenses/>.
