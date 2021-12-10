%  FILE:        end_leakage.m
%  FUNCTION:    Psl=end_leakage(rs,db,ds,dw,leo,lew,Nslts)
%  DESCRIPTION: Calculates the end leakage permeance for the inductance
%               calculations
%  INPUTS:      rs        - stator inner radius (m)
%               db        - depth of backiron (m)
%               ds        - depth of slot (m)
%               dw        - depth of winding (m)
%               leo       - length of end winding offset (m)
%               lew       - length of end winding (m)
%               Nslts     - number of slots 
%
% OUTPUTS:      Pel       - end leakage permeance (1/H)
%
% INTERNAL:     mu0       - permeability of free space (H/m)
%               rs        - stator inner radius (m)
%               wew       - width of end winding (m)   
%               Ns        - number of slots
%               lseg      - length of a segment (center of slot to center
%                           of slot) of end segment (m)
%               lr        - length of winding segment (m)
%               dr        - depth of winding segment (m)
%               lrmdr     - difference between the length and depth of
%                           winding segment (m)
%               Pel1      - permeance associated with flux flowing within 
%                           the end winding (1/H)
%               x         - temporary term associated with calculation of 
%                           Pel2 (m)
%               y         - temporary term associated with calculation of 
%                           Pel2 (m)
%               R         - reluctance of the flux path exterior to the end
%                           winding (H)
%               Pel2      - permeance associated with flux flow exterior to
%                           the end winding (1/H)
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

function Pel=end_leakage(rs,db,ds,dw,leo,lew,Nslts)

    %Handy constant
    mu0=4*pi*1e-7;
    
    % test to make sure winding present
    if (lew>0)&&(dw>0) 
     
       %Dimensions
       wew=dw;
    
       %Compute the interior permeance term for the end winding
       lseg=2*pi*(rs+ds-0.5*dw)/Nslts;
       lr=max(lew,wew);
       dr=min(lew,wew);
       lrmdr=lr-dr;
    
       if (lr==dr)
          Pel1=mu0*lseg/32; 
       else
          Pel1=mu0*lseg*(dr^4/32+ dr^3*lrmdr/16+ dr^2*lrmdr^2/64 ...
                           - dr*lrmdr^3/64+lrmdr^4*log(1+2*dr/lrmdr)/128) ...
                           /(dr^2*lr^2);
       end
    
       %Compute the exterior permeance term for the end winding
       f=@(x)((leo+lew+x)/db+ ...
              (wew+db/2+sqrt(x*(lew+leo+x)))/(2*x)+ ...
              (lew+leo+x)/(2*sqrt(x*(lew+leo+x))));
       x=fminbnd(f,1.0e-4*lew,2.0*lew);
       y=sqrt(x*(lew+leo+x));
       R=(leo+lew+x)/(db*lseg*mu0)+ ...
         (dw+0.5*db+y)/(2*x*lseg*mu0)+ ...
         (lew+leo+x)/(2*y*lseg*mu0);
       Pel2=1/R;
    
       %Compute the end leakage permeance
       Pel=2.0*(Pel1+Pel2);
       
    else
       
       Pel = 0;

    end
end