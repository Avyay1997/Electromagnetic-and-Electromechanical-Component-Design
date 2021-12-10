%  FILE:        slot_leakage.m
%  FUNCTION:    Psl=slot_leakage(g,dm,ds,l,wt,dtt,wtt,wst,wsi,dw)
%  DESCRIPTION: Calculates the slot leakage permeance for the inductance
%               calculations
%  INPUTS:      g         - air gap (m)
%               dm        - depth of permanent magnet (=0 if not surface
%                           mount PMSM machine
%               ds        - depth of slot (m)
%               l         - active length (m)
%               wt        - width of tooth (m)
%               dtt       - depth of tooth tip (m)
%               wtt       - width of tooth tip (m)
%               wst       - width of slot at tooth tip (m)
%               wsi       - width of tooth interior (m)
%               dw        - depth of winding (m)
%
% OUTPUTS:      Psl       - slot leakage permeance (1/H)
%
% INTERNAL:     pi        - 3.14159...
%               mu0       - permeability of free space (H/m)
%               dsi       - interior slot depth (m)
%               Psl1-Psl7 - permeance components associated with slot
%                           leakage permeance (1/H)     
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
% MODIFICATIONS 3/21/08. Modified by S.D. Sudhoff.  Psl1 was changed
%               based on the assumption that this flux component had
%               semicular ends.

function Psl=slot_leakage(g,dm,ds,l,wt,dtt,wtt,wst,wsi,dw)

    %Handy constant
    mu0=4*pi*1e-7;

    %Dimension
    dsi=ds-dtt;

    % Compute leakage of top of slot
    Rmax=min(dm+g,0.5*wtt);
    Psl1=mu0*l*log(1+pi*Rmax/wst)/pi;
    %  Psl1=0;   % SDS 6/16/08  Really this will interact with other slots
                 % thus am not including it
    % SDS 8/25/2010.  Decided this term really should be there          
    
    % Compute the other terms
    Psl2=l*dtt*mu0/wst;
    Psl3=l*(dsi-dw)*mu0/wsi;
    Psl4=l*dw*mu0/(3*wsi);
    Psl5=dw*mu0*log(1+pi*wt/wsi)/(3*pi);
    Psl6=mu0*log(1+pi*wt/wsi)*(dsi-dw)/pi;
    Psl7=mu0*log(1+pi*wtt/wst)*dtt/pi;
    
   
    % Find the total slot permeance
    Psl=Psl1+Psl2+Psl3+Psl4+2*(Psl5+Psl6+Psl7);
    
end