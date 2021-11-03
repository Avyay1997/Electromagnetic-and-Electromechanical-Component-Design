%  FILE:        MSE_loss.m
%  FUNCTION:    [Pld]=MSE_loss(B,qr,wr,SP)
%  DESCRIPTION: This routine finds the loss density for a steel with a 
%               given flux density waveform.  Based on a modified version
%               of the generalized Steinmetz model.

%  INPUTS:      B        - B field as electrical rotor position varied
%                          from 0 to 2*pi (T)
%               qr       - electrical rotor position.  Needs to start 
%                          at zero and go to 2*pi (rad)
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
%               N        - number of points
%               index    - vector of indices
%               pB       - the time derivative of flux density (T/s)
%               pB2      - the square of the time derivative of flux density (T/s)^2
%               Bmax     - maximum flux den+sity (T)
%               Bmin     - minimum flux density (T)
%               Bpk      - peak flux density (T)
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

function Pld=MSE_loss(B,qr,wr,SP)

   % Determine fundamental frequency
   ffund=abs(wr)/(2.0*pi);
   
   % Determine elapsed time
   t=abs(qr/wr);
   
   % Define internal parameters
   N=length(qr);
   index=[2:N-1];   
   
   % Determine time derivative of flux density
   pB(N)=(B(2)-B(N-1))/(qr(1)+2*pi-qr(N-1));
   pB(1)=pB(N);
   pB(index)=(B(index+1)-B(index-1))./(qr(index+1)-qr(index-1))';
   pB=pB*abs(wr);
   
   % Determine the square of the time derivative of flux density
   pB2=pB.^2;                     
   
   % Determine maximum B, minimum B, peak B, and range of B over a cycle
   Bmax=max(B);
   Bmin=min(B);
   Bpk=max(abs(Bmax),abs(Bmin));
   DB=Bmax-Bmin;
    
   % Determine integral
   int_dBdt2=trapz(t,pB2);
  
   % Determine equivalent frequency
   feq=2.0*int_dBdt2/(DB*pi).^2;
   
   % Determine output
   Pld=SP.MSE.kh*Bpk^(SP.MSE.beta)*feq^(SP.MSE.alpha-1)*ffund+ ...
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