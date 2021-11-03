function [rdc,rac]=resistance(ac,lc,Np,f,sigma,ur)
% resistance determines the AC and DC resistance of the wire given area,
% length, and frequency
%
% Inputs:
% ac      = cross-sectional area of copper wire (m^2)
% lc      = length of individual conductor (m)
% Np      = number of conductors in parallel (m)
% f       = vector of frequency (Hz)
% sigma   = conductor conductivity (1/(Ohm * m))
% mur     = conductor relative mermeabilty
%
% Outputs:
% rdc     = DC resistance of wire (Ohms)
% rac     = vector of AC resistancces of wire {elements corresopond to
%           f } (ohms)
%
% Internal:
% rho     = conductor resistivity (Ohm*m)
% u0      = permeability of free space (H/m)
% w       = angular frequency (rad/s)
% delta   = skin depth (m)
% r0      = fundamental resistance (Ohms)
% R       = radius of wire (m)
% R_tilda = scaled "radius" (unitless) 
% J1      = solution to 0th order Bessel function of the 1st kind (units)
% ber     = Kelvin function, real part of J1 (units)
% bei     = Kelvin function, imaginary part of J1 (units)
% deltaR  = change in scaled "radius" (unitless)
% J2      = solution to 0th order Bessel function of the 1st kind (units)
% temp    = temp. var to determine derivatives of Kelvin functions (units)
% dber    = derivative of Kelvin function, real part (units)
% dbei    = derivative of Kelvin function, imaginary part (units)

% Written by:
% Grant Shane
% Electrical Engineering Building
% 465 Northwestern Ave.
% West Lafayette, IN 47907-1285
% Phone: 765-494-3487
% E-mail: gshane@purdue.edu
% Office: EE 153

% Constants----------------------------------------------------------------
u0=pi*4e-7;   

% DC resistance of wire (Ohms)---------------------------------------------
rdc=(lc/Np)/(sigma*ac);

% AC resistance of wire (Ohms)---------------------------------------------
w=2*pi*f;                       
rho=1/sigma;   
delta=sqrt((2*rho)./(w*u0*ur));
r0=1./(sigma*delta); 
R=sqrt(ac/pi);
R_tilda=(sqrt(2)*R)./delta;   

% Bessel function
J1=besselj(0,R_tilda*exp(.75*pi*1i));

% Kelvin functions
ber=real(J1);
bei=imag(J1);

deltaR=0.0001*R_tilda;
J2=besselj(0,(R_tilda+deltaR).*exp(.75*pi*1i));
temp=(J2-J1)./deltaR;

dber=real(temp);
dbei=imag(temp);

% Resistance of wire (Ohms)
rac=(lc/Np)*(r0/(sqrt(2)*pi*R)).*(ber.*dbei-bei.*dber)./(dber.^2+dbei.^2);

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

