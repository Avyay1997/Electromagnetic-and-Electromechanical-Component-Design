%  FILE:          muB.m
%  FUNCTION:      [mu,pmu] = muB(MP,B)
%  DESCRIPTION:   This routine calculates permeability as a function of B.
%
%  Inputs:       MP       = Structure of material parameters
%                 MP.mur  = Initial relative permeability
%                 MP.muB.a[]  = Vector of alpha coefficients (1/T)
%                 MP.muB.b[]  = Vector of beta exponential coefficients (1/T)
%                 MP.muB.g[]  = Vector of gamma exponential offsets (T)
%                 MP.muB.d[]  = Vector of delta coefficients
%                 MP.muB.e[]  = Vector of epsilon values
%                 MP.muB.z[]  = Vector of zeta values
%                 MP.muB.h[]  = Vector of eta values (1/T)
%                 MP.muB.t[]  = Vector of theta values
%                B        = Vector of points at which to calculate permeability (T)
%
%  Outputs:      mu       = Permeability at points corresponding to B (H/m)
%                pmu      = derivative of mu with respect to B (H/(m*T))
%
% Internal:      aB       =  absolute value of B (T)
%                mu0      =  permeability of freespace (H/m)
%                NT       =  number of terms
%                f        =  the ratio of flux density over magnetization (with
%                            both expressd in T)
%                pg       =  derivative of g with respect to B (1/T)
%                            (g=f less a constant)
%                n        =  index variable
%                Bterm    =  term in expression
%
%  AUTHOR:      Scott D. Sudhoff                               
%               Purdue University
%               Electrical Engineering Building
%               465 Northwestern Avenue
%               West Lafayette, IN 47907-2035
%               sudhoff@ecn.purdue.edu
%               765-497-7648
%
%  REFERENCE:  Shane, G. and S. D. Sudhoff,  "Refinements in Anhysteretic
%              Characterization and Permeability Modeling,"  IEEE
%              Transactions on Magnetics, submitted for publication.

function [mu,pmu] = muB(MP,B)

%  Define some handy numbers
aB=abs(B);
mu0=4e-7*pi;
NT=length(MP.muB.a);

%  Initialize f and pg
f=MP.mur/(MP.mur-1);
pg=0;

%  Build f and pg
for n=1:NT
   Bterm=exp(-MP.muB.b(n)*aB); 
   f=f+MP.muB.a(n)*aB+MP.muB.d(n)*log(MP.muB.e(n)+MP.muB.z(n)*Bterm);
   pg=pg+MP.muB.h(n)./(MP.muB.t(n)+Bterm);
end

%  Define outputs
mu=mu0*f./(f-1);
pmu=-mu0*sign(B).*pg./(f-1).^2;

%  Copyright 2010 - Scott Sudhoff 
% 
%  The program is distributed under the terms of the GNU Lesser
%  General Public License(GNU LGPL). 
% 
%  This file is part of the Power Magnetic Material Toolbox.
% 
%  The Power Magnetic Material Toolbox is free software: 
%  you can redistribute it and/or modify
%  it under the terms of the GNU Lesser General Public License as published
%  by the Free Software Foundation, either version 3 of the License, or 
%  (at your option) any later version.
% 
%  The Power Magnetic Material Toolbox is distributed 
%  in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
%  See the GNU Lesser General Public License for more details.
% 
%  You should have received a copy of the GNU Lesser General Public
%  License along with the Power Magnetic Material Toolbox.  If not, see
%  <http://www.gnu.org/licenses/>.