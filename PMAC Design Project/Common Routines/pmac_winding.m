function [W] = pmac_winding(W,G)
% pmac_winding is a function to find the parameters related to the windings
%              of a permanent magnet ac machine.  Used in Chapter 9 of
%              "Power Magnetic Devices: A Multi-Objective Design Approach," 
%              by S.D. Sudhoff
%
% Call:
% W = pmac_winding(W,G)
%
% Inputs:
% W          = structure of winding parameters. input fields:
%  W.Ns1star = desired peak value of fundamental component 
%              of conductor density (turns/rad)
%  W.a3star  = desired ratio of third harmonic turns density to fund.                            to fundamental
%  W.apf     = packing factor (conductor area / slot area)
%  W.leo     = end winding offset (m)
% G          = structure of machine geometry (see pmac_geometry)
% 
% Outputs:
% W          = struture of winding parameters. output fields:
%  W.Ns1     = peak value of fund component of cond density (cond/rad)                            density
%  W.a3      = ratio of third harmonic turns density to fundamental
%  W.Nas     = effective series a-phase conductors in each slot (a vector)
%  W.Nbs     = effective series b-phase conductors in each slot (a vector)
%  W.Ncs     = effective series b-phase conductors in each slot (a vector)
%  W.Ns      = vector of total conductors in each slot
%  W.Mas     = effective series a-phase end conductors 
%              in each end segment (a vector)
%  W.Mbs     = effective series b-phase end conductors in 
%              each end segment (a vector)
%  W.Mcs     = effective series c-phase end conductors in  
%              each end segment (a vector)
%  W.ac      = conductor cross sectional area (m^2)
%  W.dc      = conductor diameter (m) 
%  W.dwR     = winding depth (rect. slot approx.) (m)
%  W.lew     = end winding length in direction of shaft (m)
%  W.vcd     = volume of conductor (one phase) (m^3)
%  W.vecd    = volume of end conductor (one phase) (m^3)
%  W.vscd    = volume of slot conductor (one phase) (m^3)
%             
% Internal
% r120       = 120 degrees expressed in radians
% phiss      = electrical stator position of slot centers(rad)
% hqst       = half the electrical distance of a slot and a tooth (rad)
% t1         = temporary variable for inter. calculation 
% t2         = temporary variable for inter. calculation
% t3         = temporary variable for inter. calculation
% t4         = temporary variable for inter. calculation
% Mt         = vector of total end conductors on each end segment                             segment
%
% Written by
% Brandon Cassimere and Scott D. Sudhoff                               
% Purdue University
% Electrical Engineering Building
% 465 Northwestern Avenue
% West Lafayette, IN 47907-2035
% sudhoff@purdue.edu
% 765-497-7648
               
% Physical constant
r120=2*pi/3;   

% Calculate conductor counts
phiss=G.phiss*G.P/2.0;
hqst=(pi/G.Ss)*(G.P/2);       
t1=4*W.Ns1star/G.P;
t2=sin(hqst);
t3=t1*t2;
t4=t1*W.a3star*sin(3*hqst)*sin(3*phiss)/3;
W.Nas=round(t3*sin(phiss)-t4);
W.Nbs=round(t3*sin(phiss-r120)-t4);
W.Ncs=round(t3*sin(phiss+r120)-t4);
W.Ns=abs(W.Nas)+abs(W.Nbs)+abs(W.Ncs);

% Calculation of end conductors
W.Mas=cumsum(W.Nas)-0.5*sum(W.Nas(1:(G.Ss/G.P)));                           
W.Mbs=cumsum(W.Nbs)-0.5*sum(W.Nbs(1:(G.Ss/G.P)));                           
W.Mcs=cumsum(W.Ncs)-0.5*sum(W.Ncs(1:(G.Ss/G.P))); 

% Check on end conductors
if sum(W.Mas)~=0
   error('End Winding Incorrect');
end

%Conversion from discrete description to continuous description for
%windings.  These are different from desired values because of the
%rounding which can take place
W.Ns1=W.Nas*sin(G.P*G.phiss/2)'/pi;
W.a3=-W.Nas*sin(3*G.P*G.phiss/2)'/(pi*W.Ns1);

% Compute conductor area
W.ac = G.aslt*W.apf/max(W.Ns);   % CORRECTED
W.dc = sqrt(4*W.ac/pi);        

% Compute depth of winding
W.dwR = max(W.Ns)*W.ac/(W.apf*G.wsiR);

% Compute lenght of end winding bundle
Mt=abs(W.Mas)+abs(W.Mbs)+abs(W.Mcs);
W.lew = max(Mt)*W.ac/(W.dwR*W.apf);

% Compute conductor volume
W.vscd=(G.l+2.0*W.leo)*W.ac*sum(abs(W.Nas));
W.vecd=2*pi*(G.rst+G.rsb)*W.ac*sum(abs(W.Mas))/G.Ss;
W.vcd=W.vscd+W.vecd;