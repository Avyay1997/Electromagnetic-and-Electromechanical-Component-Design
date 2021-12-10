%  FILE:        axially_compensated_length.m
%  FUNCTION:    lac=axially_compensated_length(rr,g,dm,ds,l,wt,dtt,wtt,dw,leo,lew,Ns)
%               lac=axially_compensated_length(rr,g,dm,ds,l,wt,dtt,wtt,dw,leo,lew,Ns,ri)
%  DESCRIPTION: Calculates the axially compensated length to give a more
%               accurate value for the magnetizing inductance
%  INPUTS:      rr        - rotor radius (m)
%               g         - air gap (m)
%               dm        - depth of magnet (set to zero if not surface
%                           mounted PMSM)(m)
%               ds        - depth of slot (m)
%               l         - active length (m)
%               wt        - width of tooth (m)
%               dtt       - depth of tooth tip (m)
%               wtt       - width of tooth tip (m)
%               dw        - depth of winding (m)
%               leo       - length of end winding offset (m)
%               lew       - length of end winding (m)
%               Ns        - number of slots
%               ri        - inert radius (m) {optional}
%
% OUTPUTS:      lac       - axially compensated length (m)
%
% INTERNAL:     pi        - 3.14159...
%               mu0       - permeability of free space (H/m)
%               rs        - stator inner radius (m)
%               Ns        - number of slots
%               we        - width of a stator tooth and slot segment (m)
%               de        - an approximation of the maximum distance the
%                           flux can travel before reaching the stator
%                           backiron without crossing the center of 
%                           the end winding (m)
%               w1-w4     - temporary terms associated with calculation of
%                           the permeance components asssociated with the 
%                           ends of the machine (m)
%               g1-g4     - temporary terms associated with calculation of
%                           the permeance components asssociated with the 
%                           ends of the machine (m)
%               d1-d4     - temporary terms associated with calculation of
%                           the permeance components asssociated with the 
%                           ends of the machine (m)
%               t1-t6     - temporary terms associated with calculation of 
%                           the effective length (m)
%               lseg      - length of a segment (center of slot to center
%                           of slot) of end segment (m)
%               P1-P4     - Permeance components associated with the ends
%                           of the machine (1/H)
%            
%  AUTHOR:      Brandon Cassimere                           
%               Purdue University
%               Electrical Engineering Building
%               465 Northwestern Avenue
%               West Lafayette, IN 47907-2035
%               bcassime@purdue.edu
%               765-494-3487
%
%               This code was based off of Brandon's Ph.D. 
%               dissertation "Modeling and Evolutionary Design of 
%               Permanent Magnet Sychronous Machines" (Chapter 4). 
%
% MODIFICATION  3/18/08 S.D. Sudhoff modified the calculation of de to
%               take into account conditions where either the slot depth
%               or rotor limit de.
% MODIFICATION  4/01/10 S.D.Sudhoff again modified calculation of de.
%               This involved a slightly differenct expression for the
%               nominal case; and an additional modification for the
%               case when an inert rotor is specified.

function lac=axially_compensated_length(rr,g,dm,ds,l,wt,dtt,wtt,dw,leo,lew,Ns,ri)

    %Handy constant
    mu0=4*pi*1e-7;    
    
    %Dimensions associated with the permeance calculations
    rs=rr+dm+g;
    we=2*pi*((rs+rr)/2)/Ns;    
    de1=sqrt((leo+lew/2)^2+(ds-dw/2)^2);
    de2=ds;
    if (nargin==12)
       de3=pi*rr^2/(Ns*we);
    else
       de3=pi*(rr^2-ri^2)/(Ns*we);
    end
    de=min([de1 de2 de3]);
    w1=wtt/2;
    if dtt<de
        d1=dtt;
        d3=dtt;
    else
        d1=de;
        d3=de;
    end
    g1=g+dm;
    w2=wt/2;
    d2=de-d1;
    g2=g1+2*dtt;
    w3=(we-wtt)/2;
    g3=g+dm;
    w4=(we-wt)/2;
    d4=de-d3;
    g4=g+dm+2*dtt;
        
    t1=2*pi*d3*log(1+(pi*w3)/(2*g3+2*pi*d3));
    t2=(pi*w3+2*g3)*log(1+(2*pi*d3)/(pi*w3+2*g3));
    t3=2*g3*log(1+(pi*d3)/g3);
    t4=2*pi*d4*log(1+(pi*w4)/(2*g4+2*pi*d4));
    t5=(pi*w4+2*g4)*log(1+(2*pi*d4)/(pi*w4+2*g4));
    t6=2*g4*log(1+(pi*d4)/g4);
    
    %Compute the permeances 
    P1=mu0*w1/pi*log(1+(pi*d1/g1));
    P2=mu0*w2/pi*log(1+(pi*d2/g2));
    P3=mu0/pi^2*(t1+t2-t3);
    P4=mu0/pi^2*(t4+t5-t6);

    % Compute the face permeance
    Pface=2*(P1+P2+P3+P4); 
    
    %Compute the axially compensated length
    lac=l+2*(g+dm)*Pface/(mu0*we);
    
end