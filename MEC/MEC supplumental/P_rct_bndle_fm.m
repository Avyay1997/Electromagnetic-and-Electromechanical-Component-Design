function [Pint,Pext] = P_rct_bndle_fm(ht,hb,hw,ww,cw,lw)
% P_rct_bndle_fm  Calculate the internal and external
%                 permeance associated with a rectangular
%                 bundle of conductor above a ferro or ferri - magnetic
%                 material.  Uses method similar to pg 42-43 of J. Cale,
%                 S.D. Sudhoff, L-Q Tan, "Accurately Modeling EI Core
%                 Inductors Using High-Fidelity Magentic Equivalent
%                 Circuit Approach", IEEE Transactions on Magnetics,
%                 Vol. 42, No. 1, January 2006, but accomodates clearance.
%
% [Pint,Pext] = P_rct_bndle_fm(ht,hb,hw,ww,cw,lw)
%
% Inputs:
% ht   = height of ferro/ferri magnetic material above winding (m)
% hb   = height of ferro/ferri magnetic material below winding (m)
% hw   = height of winding crossection(m)
% ww   = width of winding crossection (m)
% cw   = clearance of winding (m)
% lw   = axial length of winding bundle(m)
%
%              ->cw<-
%                 -> ww <-            
% ---------------             -
%               |             |
%               |             ht    
% magnetic      |             |
% material      |  ------     -
%               |  | w  |     |
%               |  | i  |     |
%               |  | n  |     hw
%               |  | d  |     |
%               |  | i  |     |
%               |  | n  |     |
%               |  | g  |     |
%               |  ------     -
%               |             |
%               |             hb
%               |             |
%----------------             -
%
%
% Outputs:
% P = permeance associated with the bundle (H^-1)
%
% Written by:
% S.D. Sudhoff
% School of Electrical and Computer Engineering
% 1285 Electrical Engineering Building
% West Lafayette, IN 47907-1285
% Phone: 765-494-3246
% Fax: 765-494-0676
% E-mail: sudhoff@ecn.purdue.edu

mu0=4e-7*pi;

% 1st step - calculation of reluctance of external leakage
htphb=ht+hb;
hthb=ht*hb;
Rext =((cw+ww)*htphb+sqrt(2.0*hthb*htphb*(htphb+2*hw)))/ ...
      (2.0*hthb*mu0);
Pext=1/Rext;

% 2nd step - calculation of permeance of internal leakage
k2=abs(hw-2.0*ww);
k1=k2+2.0*cw;
k3=max(ww,0.5*hw);
k4=1.25*k1-k2;
k1_2=k1*k1;
k1_3=k1_2*k1;
k3_2=k3*k3;
k3_3=k3_2*k3;
k3_4=k3_3*k3;
Pint=mu0*lw*(k3_4/3 + ...
             k4*k3_3/3 + ...
             k1*k4*k3_2/8 - ...
             k1_2*k3*k4/16 + ...
             k1_3*k4*log(1+4*k3/k1)/64)/ ...
             (ww*hw)^2;
         
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
