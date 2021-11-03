function [f]=Prob_D_EE_em_fit(P,D,fn)
% Prob_D_EE_em_fit is a fitness funtion for a design of a ui core dc biased
%           inductor, as described in Section 3.4 of "Power 
%           Magnetic Devices: A Multi-Objective Design Approach" 
%           by S.D. Sudhoff
%
% Call:
% f=Prob_D_EE_em_fit(P,D)       % 1st form - for optimization
% f=Prob_D_EE_em_fit(P,D,fn)    % 2nd form - for reporting
%
% Inputs:
% P        = parameter vector
%  P(1)    = width of base of core(m)
%  P(2)    = ratio of width of leg of core
%  P(3)    = width of conductor(m)
%  P(4)    = depth of conductor(m)
%  P(5)    = current density (A/m^2)
%  P(6)    = corner clearnace (m)
%  P(7)    = outside core clearance(m)
%  P(8)    = base core clearance

% D        = structure of design information
% D.g      = air gap (m)
% D.Bmin   = flux density (T)
% D.SZmax  = maximum allowed size sqrt(kg m^3)
% D.PEmxa  = maximum allowed loss (W)
% D.kpf    = packing factor
% D.cgd    = vertical air gap clearance(m)
% D.lc     = length of core (m)
% D.wc     = width of core (m)
% fn       = figure number (2nd form of call only)
%
% Outputs:
%  f       = vector of fitness, for 1st form of call
%  f       = structure of outputs, for 2nd form of call
%   f.SZ   = inductor size sqrt(kg m^3)
%   f.ML   = inductor mass (kg)
%   f.PL   = inductor loss (W)
%   f.J    = current density (A/m^2)
%   f.Bg   = air-gap flux density (T)
%   f.EE   = strucuture of core parameers
%   f.W    = structure of winding parameters
%   f.hL   = inductor height (m)
%   f.wL   = inductor width (m)
%   f.lL;  = inductor length (m)
%
% Internal:
% EE       = structure of EE core parameters
%  EE.MP   = structure of core material parameters
%  EE.g    = airgap (m)
%  EE.cr   = corner clearance (m)
%  EE.cso  = outside core clearance (m)
%  EE.cb   = base core clearance (m)
%  EE.csc  = clearance between center core and conductor (m)
%  EE.cgd  = air-gap clearance(m)
%  EE.lc   = length of core (m)
%  EE.we   = width of leg of E-core (m)
%  EE.wb   = width of base of E-core (m)
%  EE.wc   = width of core (m)
%  EE.ds   = slot depth (m)
%  EE.ws   = slot width (m)
%  EE.Mcr  = mass of core (kg)
% WP       = structure of wire parameters
% W        = structure of winding parameters
%  W.MP    = material parameters
%  W.ww    = width of winding bundle (m)
%  W.dw    = depth of winding bundle (m)
%  W.Vcl   = volume of coil (m^3)
%  W.Vcd   = volume of conductor (m^3)
%  W.Mcd   = mass of conductor (kg)
% ML       = inductor mass (kg)
% Bg       = air-gap flux density (T)
% PL       = inductor power loss (W)
% J        = current density (A/m^2)
% hL       = overall inductor height (m)
% wL       = overall inductor width (m)
% lL       = overall inductor length (m)
% c1       = constraint 1
% c2       = constraint 2
% c3       = constraint 3
% c4       = constraint 4
% NC       = number of constraints
% CS       = constraints satisfied
% CI       = constraints imposed 
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

% number of constraints
NC=4;

% assign parameters
J=P(5);
EE.MP=steel_catalog(5);
W.MP=conductor_catalog(1);
EE.g=D.g;
EE.cr=P(6);
EE.cso=P(7);
EE.cb=P(8);
EE.lc=D.lc;
EE.wb=P(1);
EE.we=P(2)*EE.wb;
EE.wc=D.wc;
EE.cgd=D.cgd;
EE.csc=((sqrt(2)*EE.wc+2*EE.cr)-EE.wc)/2;

% calculate winding parameters
W.ww=P(3);
W.dw=P(4);

% calculate remaining core parameters
EE.ds=W.dw+EE.cb+EE.cgd;
EE.ws=W.ww+EE.csc+EE.cso;

% compute mass
core_vol=2*(EE.we*(2*(EE.wb+EE.ds)+EE.g)*EE.lc)+...
         2*(EE.wb*(EE.wc+2*EE.ws)*EE.lc)+2*(EE.ds*EE.wc*EE.lc);
EE.Mcr=core_vol*EE.MP.row;
W.Vcl=2*pi*W.dw*(((EE.wc/2+EE.csc+W.ww)^2)-...
      ((EE.wc/2+EE.csc)^2));
W.kpf=D.kpf;
W.Vcd=W.kpf*W.Vcl;
W.Mcd=W.Vcd*W.MP.row;
ML=EE.Mcr+W.Mcd;

% compute loss
PL = J^2*W.Vcd/W.MP.sigma0;

% compute total dimensions of circumscribing box
hL=2*(EE.wb+EE.ds)+EE.g;
wL=2*(EE.we+EE.ws)+EE.wc;
lL=2*(W.ww+EE.csc)+EE.lc;
VL=hL*wL*lL;

%compute size (sqrt(kg m^3))
SZ=sqrt(ML*VL);

% check some constraints
c1=lte(SZ,D.SZmax);
c2=lte(PL,D.PLmax);

% test easy constraints
CS=c1+c2;
CI=2;
if (CS<CI)
   f=eps*[1; 1]*(CS-NC)/NC;
   return
end

% call MEC and test inductance
[Bg,c]=EE_mec(EE,W,J);
c3=mean(c);
c4=gte(Bg,D.Bmin);            

% test constraints
CS=CS+c3+c4;
CI=CI+2;
if (CS<CI)
   f=eps*[1; 1]*(CS-NC)/NC;
   return
end

% compute fitness
f=[1/SZ; 1/PL];

if (nargin>2)
   if (fn>0)
      draw_EE(EE,W,1,fn)
      disp(['EE Core Data']);
      disp(['Material = ' EE.MP.desc]);
      disp(['we (m) = ' num2str(EE.we)]);
      disp(['wc (m) = ' num2str(EE.wc)]);
      disp(['wb (m) = ' num2str(EE.wb)]);
      disp(['ws (m) = ' num2str(EE.ws)]);
      disp(['ds (m) = ' num2str(EE.ds)]);
      disp(['lc (m) = ' num2str(EE.lc)]);
      disp(['Mc (kg) = ' num2str(EE.Mcr)]);
      disp(['Winding Data']);
      disp(['Material = ' W.MP.desc]);
      disp(['dw (m) = ' num2str(W.dw)]);
      disp(['ww (m) = ' num2str(W.ww)]);
      disp(['Mcd (kg) = ' num2str(W.Mcd)]);
      disp(['Metrics']);
      disp(['SZ sqrt(kg m^3) = ' num2str(SZ)]);
      disp(['ML (kg) = ' num2str(ML)]);
      disp(['PL (W) = ' num2str(PL)]);
      disp(['J (A/m^2) = ' num2str(J)]);
      disp(['Bg (T) = ' num2str(Bg)]);
      disp(['hL (m) = ' num2str(hL)]);
      disp(['wL (m) = ' num2str(wL)]);
      disp(['lL (m) = ' num2str(lL)]);
   end
   clear f;
   f.SZ=SZ;
   f.ML=ML;
   f.VL=VL;
   f.PL=PL;
   f.J=J;
   f.Bg=Bg;
   f.EE=EE;
   f.W=W;
   f.hL=hL;
   f.wL=wL;
   f.lL=lL;
end % if (nargin>2)
    
end

function c=lte(x,xmx)
% LTE      less-than-or-equal to function
%
% c=lte(x,xmx)
%
% Inputs:
% x        = a quantitity 
% xmx      = maximum allowed value
%
% Outputs:
% c        = constraint variable - 1 if x<=xmx, 0<c<1 if x>xmx

   if (x<=xmx)
      c=1;
   else
      c=1/(1+x-xmx);
   end
   
end

function c=gte(x,xmn)
% GTE        greater-than-or-equal to function
%
% c=gte(x,xmn)
%
% Inputs:
% x        = a quantitity 
% xmn      = minimum allowed value
%
% Outputs:
% c        = constraint variable - 1 if x>=xmn, 0<c<1 if x<xmn

   if (x>=xmn)
      c=1;
   else
      c=1/(1+xmn-x);
   end

end