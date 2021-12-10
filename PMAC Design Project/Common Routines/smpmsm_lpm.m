%  FILE:        smpmsm_lpm.m
%  FUNCTION:    [Rs,Lq,Ld,lpm]=smpmsm_lpm(P,rr,g,dm,db,ds,l, ...
%                                         wt,dtt,wtt,wst,wsi,dw, ...
%                                         sc,ac,leo,lew, ...
%                                         apm,Bpm,murm,murs, ...
%                                         Np,Nas,Mas,Nbs,Mbs)
%               [Rs,Rss,Rse,vct,vce,vce,Lq,Ld,lpm]=smpmsm_lpm( ...
%                                         P,rr,g,dm,db,ds,l, ...
%                                         wt,dtt,wtt,wst,wsi,dw, ...
%                                         sc,ac,leo,lew, ...
%                                         apm,Bpm,murm,murs, ...
%                                         Np,Nas,Mas,Nbs,Mbs)
%               [Rs,Lq,Ld,lpm]=smpmsm_lpm(P,rr,g,dm,db,ds,l, ...
%                                         wt,dtt,wtt,wst,wsi,dw, ...
%                                         sc,ac,leo,lew, ...
%                                         apm,Bpm,murm,murs, ...
%                                         Np,Nas,Mas,Nbs,Mbs,ri)
%               [Rs,Rss,Rse,vct,vce,vce,Lq,Ld,lpm]=smpmsm_lpm( ...
%                                         P,rr,g,dm,db,ds,l, ...
%                                         wt,dtt,wtt,wst,wsi,dw, ...
%                                         sc,ac,leo,lew, ...
%                                         apm,Bpm,murm,murs, ...
%                                         Np,Nas,Mas,Nbs,Mbs,ri)
%  DESCRIPTION: Calculates lumped parameter model parameters of a
%               permanent magnet synchronous machine
%  INPUTS:      P         - number of poles (should be even)
%               rr        - rotor radius (m)
%               g         - air gap (m)
%               dm        - depth of magnet (m)
%               db        - depth of backiron (m)
%               ds        - depth of slot (m)
%               l         - active length (m)
%               wt        - width of tooth (m)
%               dtt       - depth of tooth tip (m)
%               wtt       - width of tooth tip (m)
%               wst       - width of slot at tooth tip (m)
%               wsi       - width of slot interior (m)
%               dw        - depth of winding (m)
%               sc        - conductivity of conductors (mhos/m^2)
%               ac        - effective series conductor cross sectional area
%                           (m^2);
%               leo       - length of end winding offset (m)
%               lew       - length of end winding (m)
%               apm       - permanet magnet fraction (m)
%               Bpm       - residual flux density (T)
%               murm      - relative permeability of permanent magnet 
%               murs      - relative permeability of spacer
%               Np        - peak value of fund. component of turns density
%               Nas       - eff. series a-phs cond. in each slot (a vector)
%               Nbs       - eff. series b-phs cond. in each slot (a vestor)
%               Mas       - eff. series a-phs end conductors in each end 
%                           segment (a vector)
%               Mbs       - eff, series b-phase end conductors in each end 
%                           segment (a vector)
%               ri        - inert radius (m) (optional)
%
% OUTPUTS:      Rs        - stator winding resistance (Ohms)
%               Rss       - stator wind. resist. assoc. with slots (Ohms)
%               Rse       - stat. wind. resist. assoc. with both end (Ohms)
%               vct       - total conductor volume (m^3)
%               vcs       - slot conductor volume (m^3)
%               vce       - slot end volume (m^3)
%               Lq        - q-axis inductance (H)
%               Ld        - d-axis inductance (H)
%               lm        - flux linkage due to permanent magnet (Vs)
%
% INTERNAL:     pi        - 3.14159...
%               mu0       - permeability of free space (H/m)
%               geff      - effective air gap (m)
%               rs        - stator inner radius (Ohms)
%               Nslts     - number of slots
%               Rsp       - quasi-reluctance associated with spacer (1/H)
%               Rpm       - quasi-reluctance assoc. with permanent magnet
%                           (1/H)
%               Rg        - quasi-reluctance associated with air gap (1/H)
%               lac       - effective length for use in calcuation of
%                           magnetizing inductances (H)
%               Lqm       - q-axis magnetizing inductance (H)
%               Ldm       - d-axis magnetizing inductance (H)
%               Psl       - Slot leakage permeance (1/H)
%               Pel       - end leakage permeance (1/H)
%               Llp       - self leakage inductance (H)
%               Llm       - mutual leakage inductance (H)
%               Lls       - stator leakage inductance (H)
%
% Modifications
%   Modified on 2/11/2010 by S.D. Sudhoff.  Code adjusted to provide
%   breakout of end turn and slot stator resistance and respective
%   volumes
    
function [varargout]=smpmsm_lpm(P,rr,g,dm,db,ds,l, ...
                                wt,dtt,wtt,wst,wsi,dw, ...
                                sc,ac,leo,lew, ...
                                apm,Bpm,murm,murs, ...
                                Np,Nas,Mas,Nbs,Mbs,ri)
    % handy constants
    mu0=4*pi*1e-7;
    
    % compute the stator radius
    rs=rr+g+dm;
    
    % compute the number of slots
    Nslts=length(Nas);
    
    %Compute the effective air gap using Carter's coefficient
    geff=effective_airgap(g,wt,dtt,wtt,wst);
    
    % compute quasi-reluctances
    Rsp=rr*log(1+dm/rr)/(mu0*murs);
    Rpm=rr*log(1+dm/rr)/(mu0*murm);
    Rg =rr*log(1+geff/(rr+dm))/mu0;
 
    % compute flux density due to the permanent magnet
    lm=8.0*rr*l*Np*dm*Bpm*sin(0.5*pi*apm)/(P*(Rpm+Rg)*mu0*murm);
    
    % Compute the axially compensated length
    if (nargin==26)
       lac=axially_compensated_length(rr,g,dm,ds,l,wt,dtt,wtt,dw, ...
                                      leo,lew,Nslts);
    else
       lac=axially_compensated_length(rr,g,dm,ds,l,wt,dtt,wtt,dw, ...
                                      leo,lew,Nslts,ri);
    end
    
    % compute the magnetizing inductances    
    Lqm=6*lac*rr*Np^2*( (pi*(1-apm)+sin(pi*apm))/(Rsp+Rg) + ...
                         (pi*apm -   sin(pi*apm))/(Rpm+Rg))/P^2;
    Ldm=6*lac*rr*Np^2*( (pi*(1-apm)-sin(pi*apm))/(Rsp+Rg) + ...
                         (pi*apm +   sin(pi*apm))/(Rpm+Rg))/P^2;
    
    % Compute slot leakage permeance
    Psl=slot_leakage(g,dm,ds,l,wt,dtt,wtt,wst,wsi,dw);
        
    % compute end leakage permeance
    Pel=end_leakage(rs,db,ds,dw,leo,lew,Nslts);
           
    % compute leakage inductances
    Llp=Psl*sum(Nas.^2)+Pel*sum(Mas.^2);
    Llm=Psl*sum(Nas.*Nbs)+Pel*sum(Mas.*Mbs);
    Lls=Llp-Llm;
    
    % compute the q- and d-axis inductances
    Lq=Lqm+Lls;
    Ld=Ldm+Lls;
    
    % compute the stator resistance
    [Rs,Rss,Rse,vct,vcs,vce]=stator_resistance(rs,ds,dw,l,leo,lew, ...
                                               ac,sc,Nas,Mas);
    
    % assign output arguments
    if nargout==4
       varargout={Rs,Lq,Ld,lm};
    else
       varargout={Rs,Rss,Rse,vct,vcs,vce,Lq,Ld,lm};
    end

    if (Lq<0)
        disp('Debugging smpmsm_lpm');
        keyboard
    end
end