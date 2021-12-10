%  FILE:        stator_resistance.m
%  FUNCTION:    [Rs]=stator_resistance(rs,ds,dw,l,leo,lew,ac,sc,Nas,Mas)
%               [Rs,Rss,Rse,vct,vcs,vce] = 
%                  stator_resistance(rs,ds,dw,l,leo,lew,ac,sc,Nas,Mas)
%  DESCRIPTION: Calculates the stator resistance
%  INPUTS:      rs        - stator innner radius (m)
%               ds        - depth of slot (m)
%               dw        - depth of winding (m)
%               l         - active length (m)
%               leo       - length of end winding offset (m)
%               lew       - length of end winding (m)
%               sc        - conductivity of conductors (mhos/m^2)
%               ac        - effective series conductor cross sectional area
%                           (m^2)
%               Nas       - effective series a-phase conductors in each
%                           slot (a vector)
%               Mas       - effective series a-phase end conductors 
%                           in each end segment (a vector)
%
% OUTPUTS:      Rs        - stator resistance (Ohms)
%               Rss       - stator resistance contributions of slots (Ohms)
%               Rse       - stator resistance contributions of ends (Ohms)
%                           note: this includes both ends
%               vct       - total a-phase conductor volume (m^3)
%               vcs       - volume of a-phase conductor in slots (m^3)
%               vce       - volume of a-phase conductor in end turns (m^3)
%                           note: this includes both ends
%
% INTERNAL:     Ns        - number of slots
%               lseg      - length of segment of end turns (m)
%               ac2s      - conductor area ^2 * conductivity (m/Ohm)
%               Nsc       - number of slot conductors per phase
%
%  AUTHOR:      Brandon Cassimere                           
%               Purdue University
%               Electrical Engineering Building
%               465 Northwestern Avenue
%               West Lafayette, IN 47907-2035
%               bcassime@purdue.edu
%               765-494-3487
%
%               This code was based off of a monograph by S.D. Sudhoff
% Modifications
%   Modified on 2/11/2010 by SDS to seperate end turn and slot resistance

function [varargout]=stator_resistance(rs,ds,dw,l,leo,lew,ac,sc,Nas,Mas)

    %Dimensions
    Ns=length(Nas);
    lseg=2*pi*(rs+ds-0.5*dw)/Ns;
    Nsc=sum(abs(Nas));
    
    %Compute volumes
    vcs =l*ac*Nsc;
    vce =2*ac*(lseg*sum(abs(Mas))+(leo+0.5*lew)*Nsc);
    vct =vcs+vce;
    
    % Compute resistances
    ac2s=sc*ac^2;
    Rss=vcs/ac2s;
    Rse=vce/ac2s;
    Rs=Rss+Rse;
    
    if nargout==1
       varargout={Rs};
    else
       varargout={Rs,Rss,Rse,vct,vcs,vce};
    end
    
end