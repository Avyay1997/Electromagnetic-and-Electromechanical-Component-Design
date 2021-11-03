function [mu,pmu] = muB(MP,B)

% muB        Calculates permeability as a function of B (flux density)
%
% [mu,pmu]  = muB(MP,B)
%
% Inputs:
% MP       = Structure of material coefficeints
%  MP.mur  = Initial relative permeability
%  MP.muB.a[]  = Vector of alpha coefficients (H/(m*T))
%  MP.muB.b[]  = Vector of beta coefficients (1/T)
%  MP.muB.g[]  = Vector of gamma coefficients (T)
%  MP.muB.d[]  = Vector of delta coefficients (H/m)
%  MP.muB.e[]  = Vector of epsilon values
%  MP.muB.z[]  = Vector of zeta values
%  MP.muB.h[]  = Vector of eta values (H/(m*T))
%  MP.muB.t[]  = Vector of theta values
% B        = Vector of point at which to calculate permeability
%
% Outputs:
% mu       = Permeability at points corresponding to B (H/m)
% pmu      = derivative of mu with respect to B (H/(m*T))
%
% Internal:
% aB       =  absolute value of B
% f        =  the ratio of flux density over magnetization (with
%             both expressd in T)
% g        =  f less a constant
% pg       =  derivative of g with respect to B
% n        =  index variable
% NT       =  number of terms
% mu0      =  permeability of freespace
%
% Written by:
% S.D. Sudhoff
% School of Electrical and Computer Engineering
% 1285 Electrical Engineering Building
% West Lafayette, IN 47907-1285
% Phone: 765-494-3246
% Fax: 765-494-0676
% E-mail: sudhoff@ecn.purdue.edu


aB=abs(B);
mu0=4e-7*pi;
NT=length(MP.muB.a);
f=MP.mur/(MP.mur-1);
pg=0;
for n=1:NT
   Bterm=exp(-MP.muB.b(n)*aB); 
   f=f+MP.muB.a(n)*aB+MP.muB.d(n)*log(MP.muB.e(n)+MP.muB.z(n)*Bterm);
   pg=pg+MP.muB.h(n)./(MP.muB.t(n)+Bterm);
end
mu=mu0*f./(f-1);
pmu=-mu0*sign(B).*pg./(f-1).^2;


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
%  License along with MEC Toolbox.  
%  If not, see <http://www.gnu.org/licenses/>.

