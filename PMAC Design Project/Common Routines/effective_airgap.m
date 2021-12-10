%  FILE:        effective_airgap.m
%  FUNCTION:    geff=effective_airgap(g,wt,dtt,wtt,wst)
%  DESCRIPTION: Calculates the effective air gap using Carter's coefficient
%  INPUTS:      g         - air gap (m)
%               wt        - width of tooth (m)
%               dtt       - depth of tooth tip (m)
%               wtt       - width of tooth tip (m)
%               wst       - width of slot at tooth tip (m)
%
% OUTPUTS:      geff      - effective air gap (m)

% INTERNAL:     pi        - 3.14159...
%               ralpha    - minimum distance between the depth of tooth tip
%                           and half the width of slot at tooth tip (m)
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

function geff=effective_airgap(g,wt,dtt,wtt,wst)

% Compute the effective air using Carter's coefficient
    ralpha=min(dtt,0.5*wst);
    geff=g*(wst+wtt)/(wtt+4*g*(log(1+0.5*pi*ralpha/g)+ ...
          log((4*g+2*(wtt-wt)+pi*wst)/ ...
          (4*g+2*(wtt-wt)+2*pi*ralpha)))/pi);

end
