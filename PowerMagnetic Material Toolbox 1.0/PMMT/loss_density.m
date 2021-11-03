%  FILE:        loss_density.m
%  FUNCTION:    [Pld]=loss_density(B,pB,wr,SP)
%  DESCRIPTION: This routine finds the loss density for a steel with a 
%               given flux density waveform.  Based on a modified version
%               of the generalized Steinmetz model.  NOTE THAT THE INPUT
%               WAVEFORMS NEED TO BE OVER AN INTEGER # CYCLES.

%  INPUTS:      B        - B field at corresponding electrical
%                          rotor positions (data over int. # cycles) (T)
%               pB       - derivative of B field at corresponding
%                          rotor positions (with respect to electrical
%                          rotor position) (data over int. # cycles)(T/rad)
%               wr       - electrical rotor speed (rad/s)
%               SP       - structure of steel parameters
%                 SP.MSE.alpha  - frequency exponent
%                 SP.MSE.beta   - flux density exponent
%                 SP.MSE.kh     - hysteresis loss coefficient (J/m^3)
%                 SP.MSE.ke     - eddy current loss coefficient (J*s/m^3)
%
%  OUTPUTS:     Pld      - power loss density (W/m^3)
%
%  INTERNAL:    t        - time since start of a cycle (s)
%               ffund    - fundamental frequency (Hz)
%               pB2      - the square of the time derivative of flux density (T/s)^2
%               Bmax     - maximum flux density (T)
%               Bmin     - minimum flux density (T)
%               DB       - range of flux density over a cycle (T)
%               int_dBdt2- integral of square of time derivative of flux density
%                          over a cycle (T^2/s)
%               feq      - equivalent frequency (Hz)
%
%  AUTHOR:      Scott D. Sudhoff                               
%               Purdue University
%               Electrical Engineering Building
%               465 Northwestern Avenue
%               West Lafayette, IN 47907-2035
%               sudhoff@ecn.purdue.edu
%               765-497-7648
%
%  REFERENCE:   Sudhoff, S. D., ECE 695 (Lecture Set 4: Modeling Magnetic
%               Materials).  West Lafayette, IN: Purdue University, Dept. 
%               of Electrical and Computer Engineering, 2009.

function Pld=loss_density(B,pB,wr,SP)

   % Determine fundamental frequency
   ffund=wr/(2.0*pi); 
   
   % Define range of time
   t=linspace(0,1/ffund,length(B));
   
   % Determine the square of the time derivative of B
   pB2=abs(pB*wr).^2;                     
   
   % Determine maximum B, minimum B, and range of B over a cycle
   Bmax=max(B);
   Bmin=min(B);
   DB=Bmax-Bmin;
    
   % Determine integral
   int_dBdt2=trapz(t,pB2);
  
   % Determine equivalent frequency
   feq=2.0*int_dBdt2/(DB*pi).^2;
   
   % Determine output
   Pld=SP.MSE.kh*Bmax^(SP.MSE.beta)*feq^(SP.MSE.alpha-1)*ffund+ ...
       SP.MSE.ke*ffund*int_dBdt2;
   
end

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