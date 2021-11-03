%  FILE:          resistance.m
%  FUNCTION:      [rdc,rac]=resistance(ac,lc,Np,f,sigma,ur)
%  DESCRIPTION:   This routine determines the AC and DC resistances of a 
%                 wire.  Note the resistances are calculated
%                 based on an insulated, solid conductor within the 
%                 wire.  A wire may consist of multiple insulated, solid 
%                 conductors.
%
%  INPUTS:      ac      = cross-sectional area of an insulated, solid
%                         conductor within the wire (m^2)
%               lc      = length of conductor (m)
%               Np      = number of insulated, solid conductors in parallel
%                         that form the wire
%               f       = vector of frequencies (Hz)
%               sigma   = conductivity of conductor (1/(Ohm*m))
%               mur     = relative permeabilty of conductor
%
%  OUTPUTS:
%               rdc     = DC resistance of wire (Ohms)
%               rac     = vector of AC resistances of wire (ohms)
%                         {elements correspond to f} 
%
%  INTERNAL:
%               rho     = conductor resistivity (Ohm*m)
%               u0      = permeability of free space (H/m)
%               w       = angular frequency (rad/s)
%               delta   = skin depth (m)
%               r0      = fundamental resistance (Ohms)
%               R       = radius of wire (m)
%               R_tilda = scaled "radius" (unitless) 
%               J1      = solution to 0th order Bessel function of the 1st 
%                         kind (units)
%               ber     = Kelvin function, real part of J1 (units)
%               bei     = Kelvin function, imaginary part of J1 (units)
%               deltaR  = change in scaled "radius" (unitless)
%               J2      = solution to 0th order Bessel function of the 1st 
%                         kind (units)
%               temp    = temp. var to determine derivatives of Kelvin 
%                         functions (units)
%               dber    = derivative of Kelvin function, real part (units)
%               dbei    = derivative of Kelvin function, imaginary part (units)
%
%  AUTHOR:      Grant Shane for Scott D. Sudhoff                               
%               Purdue University
%               Electrical Engineering Building
%               465 Northwestern Avenue
%               West Lafayette, IN 47907-2035
%               sudhoff@ecn.purdue.edu
%               765-497-7648
%
%  REFERENCE:   Ramo, S., J. R. Whinnery, and T. Van Duzer, Fields and 
%               Waves in Communication Electronics.  New York: John Wiley &
%               Sons, Inc., 1965.

function [rdc,rac]=resistance(ac,lc,Np,f,sigma,ur)

%  Constant
u0=pi*4e-7;   

%  Determine DC resistance of wire
rdc=(lc/Np)/(sigma*ac);

%  Parameters needed to determine AC resistance of wire
w=2*pi*f;                       
rho=1/sigma;   
delta=sqrt((2*rho)./(w*u0*ur));
r0=1./(sigma*delta); 
R=sqrt(ac/pi);
R_tilda=(sqrt(2)*R)./delta;   

%  Bessel function
J1=besselj(0,R_tilda*exp(.75*pi*1i));

%  Kelvin functions
ber=real(J1);
bei=imag(J1);

%  Define change in R_tilda
deltaR=0.0001*R_tilda;

%  Bessel Function
J2=besselj(0,(R_tilda+deltaR).*exp(.75*pi*1i));

%  Determine slope
temp=(J2-J1)./deltaR;

%  Derivatives of Kelvin functions
dber=real(temp);
dbei=imag(temp);

%  AC Resistance of wire
rac=(lc/Np)*(r0/(sqrt(2)*pi*R)).*(ber.*dbei-bei.*dber)./(dber.^2+dbei.^2);

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