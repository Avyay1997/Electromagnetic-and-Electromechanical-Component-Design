%  FILE:        smpmsm_field_analysis3.m
%  FUNCTION:    [Phit,Phib,Bt,Bb,Btb,qr,Btqr,Bbqr, ...
%               Bt_max,Bb_max,Btb_max,Brot_max,Hm_min,Ptcl,Pbcl]= ...
%               smpmsm_field_analysis3(Nrp,P,ns,l,rr,rm,rb,ro,wt,db, ...
%                                        rf,apt,np,dm,alphapm,bpm, ...
%                                        mur_pm,mur_sp,g,is,phii,wr,SS,fn)
%  DESCRIPTION: permanent magnet synchronous machine analysis
%  INPUTS:      Nrp      - number of rotor positions
%               P        - number of poles
%               ns       - number of slots
%               l        - active length (m)
%               rr       - rotor iron radius (m)
%               ri       - rotor iron inert radius (m)
%               rsh      - shaft radius (m)
%               rm       - radius to outside edge of PM (m)
%               rb       - radius to inside edge of backiron (m)
%               ro       - radius to outside edge of backiron (m)
%               wt       - width of tooth (m)
%               db       - depth of backiron (m)
%               rf       - fillet radius (m)
%               apt      - principal tooth area (m^2)
%               np       - actual peak value of fundamental component  
%                             of conductor density (turns/rad)
%               dm       - magnet depth (m)
%               alphapm  - permanent magnet fraction
%               bpm      - residual magnetization of PM (T)
%               mur_pm   - relative permeability of PM 
%               mur_sp   - relative permeability of spacer
%               g        - air gap (m)
%               is       - RMS phase current (A)
%               phii     - current phase advance (rad)
%               wr       - electrical rotor speed (rad/s)
%               SS       - steel desription 
%               fn       - function handle to output the report
%
% OUTPUTS:      Phit     - number of rotor positions by number of slots
%                          vector of flux in each tooth for each rotor 
%                          position (Wb)
%               Phib     - number of rotor positions by number of slots
%                          vector of flux in each backiron segments (Wb)
%               Bt       - number of rotor positions by number of slots
%                          vector of flux density in each tooth (T)
%               Bb       - number of rotor positions by number of slots
%                          vector of flux density in each segment (T)
%               Btb      - number of rotor positions by number of slots
%                          vector of flux density at the tooth base (T)
%               Bt_max   - maximum tooth flux density (T)
%               Bb_max   - maximum backiron flux density (T)
%               Btb_max  - maximum tooth base flux density (T)
%               Brot_max - maximum flux density on the rotor (T) 
%               Hm_min   - minimum field intensity in region of positive
%                          magnitization.  This may be used to test
%                          whether or not demagnetization could occur
%                          (A/m)
%               Ptcl     - tooth core (hysteresis and eddy) loss (W)
%               Pbcl     - backiron core (hysteresis and eddy) loss (W)
%
%  INTERNAL:    mu0         - permeability of free space (H/m)
%               ws          - width of slot (m)
%               cs          - Carter's coeffiecient
%               ge          - effective air gap (m)
%               r120        - 120 degrees expressed in radians
%               pi          - 3.14159....
%               qrm         - mechanical rotor position considered in 
%                             analysis (rad)
%               Fpm         - MMF due to PM in PM region (A)
%               Rpm         - reluctance density of PM region ((H/m^2)^-1)
%               Rsp         - reluctance density of spacer reg ((H/m^2)^-1)
%               Rg          - reluctance density of air gap reg((H/m^2)^-1)
%               rp          - rotor position index
%               t1-t2       - terms for calculations
%               phirm_mt    - mechanical position of PM transition
%                             points relative to rotor (rad)
%               phism_mt    - electrical position of PM transition points
%                             relative to stator.  Transistion points are
%                             those point where we have a change in    
%                             magnetization of the PM; or on effective
%                             tooth boundaries (rad)
%               phism_tp    - electrical position of all transition points
%                             relative to stator (PM + tooth edge) (rad)
%               region_array- an array whose elements indicate the PM reg.
%                             of each interval between the trans. points.
%                             The PM regions are as follows.  Region 1:
%                             positive magnetization.  Region 2: spacer.
%                             Region 3: Negative magnetization.  Region 4:
%                             spacer.  Region 5: Positive magnetization.
%               tooth_array - an array whose elements indicate which tooth
%                             each interval between trans. points is in.
%               tooth       - tooth number
%               region      - PM region number
%               valid_list  - list of intervals in first poles worth of 
%                             teeth
%               Fm          - Array whose elements are the MMF due to the
%                             PM in each interval (A)
%               Rp          - Array whose elements are the reluctance 
%                             density ((H/m^2)^-1)
%               pos_mag     - Array of index where the positive 
%                             magnetization occurs
%               neg_mag     - Array of index where the negative 
%                             magnetization occurs
%               phism_sc    - mechanical stator position of slot centers 
%                             (rad)
%               phism_tc    - mechanical stator position of tooth centers 
%                             (rad)
%               phism_te    - mechanical stator position of tooth edge
%                             (leading edge, were we start counting flux,
%                              not the actual location of the tooth) (rad)
%               index1      - index of beginning of integration interval
%               index2      - index of end of integration interval
%               phism1      - mechanical stator position at beginning of
%                             integration interval (rad)
%               phism2      - mechanical stator position at end of 
%                             integration interval (rad)
%               phir1       - electrical rotor position at beginning of
%                             integration interval (rad)
%               phir2       - electrical rotor position at end of 
%                             integration interval (rad)
%               Fs1         - Stator MMF at the beginning of integration
%                             (A)
%               Fs2         - Stator MMF at the end of integration (A)
%               phir        - vector of electrical rotor position points
%                             (rad)
%               phir_min    - electrical rotor position where stator MMF is
%                             minimum (rad)
%               phir_max    - electrical rotor position where stator MMF is
%                             maximum (rad)
%               Fs_min      - minimum value of the stator MMF (A)
%               B1-B5       - rotor flux density at candidate positions (T)
%               qr          - electrical rotor positions (rad)
%               pBtqr       - derivative of B field on the tooth with 
%                             respect to electrical rotor position 
%                             at corresponding rotor positions (T/rad) 
%               Btqr        - B field on the tooth at corresponding 
%                             electrical rotor positions (T)
%               Ptcl        - Tooth core losses (hysteresis and 
%                             eddy currents) (W)
%               pBbqr       - derivative of B field on the back-iron with 
%                             respect to electrical rotor position 
%                             at corresponding rotor positions (T/rad) 
%               Bbqr        - B field on the back-iron at corresponding 
%                             electrical rotor positions (T)
%               Pbcl        - Back-iron core losses (hysteresis and 
%                             eddy currents) (W)
%               Brot_max1   - maximum rotor flux density based on
%                             transition points (T)
%               Brot_rad    - maxtrix of number of rotor positions by 
%                             number of slots with average radial flux
%                             desnsity in rotor over tooth/slot region (T)
%               Brot_tan    - maxtrix of number of rotor positions by 
%                             number of slots with average tangential flux
%                             desnsity in rotor over tooth/slot region (T)
%               Brot_mag    - maxtrix of number of rotor positions by 
%                             number of slots with average magnitude flux
%                             desnsity in rotor over tooth/slot region (T)
%               Brot_max2   - max rotor flux density based on Brot_mag (T)
%
%  AUTHOR:      Scott D. Sudhoff                               
%               Purdue University
%               Electrical Engineering Building
%               465 Northwestern Avenue
%               West Lafayette, IN 47907-2035
%               sudhoff@ecn.purdue.edu
%               765-497-7648
%
%               Parts of this code were based off of an earlier version
%               (pmsm_analysis) written by J. Cale, based on a monograph
%               manuscript by S. Sudhoff.  Brandon Cassimere had an
%               important role in both this code and (pmsm_analyis) in
%               identifying and correcting numerous errors.
%
%               This code is different from smpmsm_field_analysis2 in that
%               it includes Carter's coefficient

function [Phit,Phib,Bt,Bb,Btb,qr,Btqr,Bbqr,Bt_max,Bb_max,Btb_max, ...
          Brot_max,Hm_min,Ptcl,Pbcl]= ...
          smpmsm_field_analysis3(Nrp,P,ns,l,rr,ri,rsh,rm,rb,ro,wt,db, ...
                                rf,apt,np,dm,alphapm,bpm, ...
                                mur_pm,mur_sp,g,is,phii,wr,SS,phism_te,phism_tc,fn)

% Physical constants
mu0=4*pi*10^-7;         
r120=2*pi/3;            


% compute effective air gap
ws=2*pi*(rr+dm+g)/ns-wt;
cs=(wt+ws)/(wt+4*g*log(1+0.25*pi*ws/g)/pi);
ge=g*cs;

% Field Solution----------------------------------------------------------%

% Rotor positions considered
if Nrp>1
  qrm=(pi/ns)*linspace(-1,(Nrp-2)/Nrp,Nrp);            
else
  qrm=0.0;
end 
thetarm=qrm;

% Compute magnetizing inductances and flux linkage due to PM
Fpm=dm*bpm/(mu0*mur_pm);
Rpm=rr*log(1+dm/rr)/(mur_pm*mu0);
Rsp=rr*log(1+dm/rr)/(mur_sp*mu0);
Rg =rr*log(1+ge/(rr+dm))/mu0;

% The following code computes the tooth flux density
% The operation of the code is as follows.  First
% a set of PM transition points are defined (phirm_mt).
% These are the points where the PM has a change in magnetization.
% There are also transition points on the effective tooth edge.
% The effective tooth edge denote the boundaries of integration
% when determining the flux in the tooth.
% Next, all transition points are stacked together and sorted, and
% the intervals between transition points are identified.
% Each of these regions is within a single tooth flux collection
% area and has a constant magnetization.
% The list of intervals is then truncated to correspond to one
% pole's worth of teeth.  
% Next, the flux in each of these teeth is determined, and the
% remaining teeth are found by symmetry.

   % Compute PM transition locations
   t1=pi*(1-alphapm)/2;
   t2=2/P;
   phirm_mt=[-t1 t1 pi-t1 pi+t1 2*pi-t1]*t2;

   % Initialize flux density in tooth
   Phit=zeros(Nrp,ns);

   % Calculate the flux in the first pole's worth of teeth
   for rp=1:Nrp
    
      % compute transition angles for PM material in terms of stator
      % position
      phism_mt=phirm_mt+qrm(rp);
    
      % make list of transition points (either because tooth changes,
      % or because PM changes
      phism_tp=sort([phism_mt phism_te]);
   
      % create an array describing region and tooth
      region_array=ones(1,length(phism_tp));
      
      tooth_array =ones(1,length(phism_tp));
      tooth=0;
      region=1;
      for tp=1:length(phism_tp)
         tooth_index=tooth+1;
         if (tooth_index>ns)
            tooth_index=1;
         end
         if phism_tp(tp)>=phism_te(tooth_index)
            tooth=tooth+1;
         end
         if (phism_tp(tp)>=phism_mt(region))&(region<=4)
            region=region+1;
         end
         region_array(tp)=region;
         tooth_array(tp)=tooth;
      end
             
      % truncate list to one pole's worth of teeth
      valid_list   = find((phism_tp>=phism_te(1))& ...
                          (phism_tp<=phism_te(1+ns/P)));
      region_array = region_array(valid_list);
      tooth_array  = tooth_array(valid_list);
      phism_tp     = phism_tp(valid_list);
    
      % assign magnetization and reluctance between each point
      Fm=zeros(1,length(phism_tp)-1);
      Rp=ones(1,length(phism_tp)-1)*Rsp;
      pos_mag=find((region_array==1)|(region_array==5));
      neg_mag=find(region_array==3);
      Fm(pos_mag)= Fpm;
      Fm(neg_mag)=-Fpm;
      Rp(pos_mag)=Rpm;
      Rp(neg_mag)=Rpm;
      
      % integrate away
      t1=6.0*sqrt(2.0)*np*is/P^2;
      for j=2:length(phism_tp)
         index1=j-1;
         index2=j;
         phism1=phism_tp(index1);
         phism2=phism_tp(index2);
         phir1=P*(phism1-qrm(rp))/2.0;
         phir2=P*(phism2-qrm(rp))/2.0;
         tooth=tooth_array(index1);
         Phit(rp,tooth)=Phit(rp,tooth)+ ...
                          (rr*l/(Rp(index1)+Rg))* ...
                          (Fm(index1)*(phism2-phism1)+ ...
                           t1*(sin(phir2-phii)- ...
                               sin(phir1-phii)));
      end
    
      % compute the flux in the rest of the teeth
      t1=ns/P;
      for j=t1+1:ns
         Phit(rp,j)=-Phit(rp,j-t1);
      end
    
   end
   
% End calculation of the tooth flux

% Compute the backiron flux
Phib(Nrp,ns)=0;
Phib(:,1)=Phit(:,1)-0.5*sum(Phit(:,1:ns/P),2);
for j=2:ns
    Phib(:,j)=Phib(:,j-1)+Phit(:,j);
end

% Compute the statpr backiron and tooth flux density
Bt=Phit/(l*wt);
Bb=Phib/(l*db);

% Calculation of the flux density at the base of the ith tooth
t1=wt/(wt+2.0*rf);
Btb(Nrp,ns)=0;
Btb(:,1)=((Bt(:,1)*t1).^2+0.25*(Bb(:,1)+Bb(:,ns)).^2).^0.5;
for i=2:ns 
   Btb(:,i)=((Bt(:,i)*t1).^2+0.25*(Bb(:,i)+Bb(:,i-1)).^2).^0.5;
end


%Find the minimum H in region of positive magnetization (by symmetry, then
%we don't need to consider the region of negative magnetization
%Note, the minimum point must either be at the transition points, or the
%point where the stator MMF is a minimum (if this point falls in the pos.
%PM magnetization region
t1=3*sqrt(2.0)*np*is/P;
t2=(1-alphapm)*pi/2;
Fs1=t1*cos(pi+t2-phii);
Fs2=t1*cos(2*pi-t2-phii);
phir_min=pi+phii;
if (phir_min>2*pi)
   phir_min=phir_min-2*pi;
end
if (phir_min<0)
   phir_min=phir_min+2*pi;
end
if ((phir_min>pi+t2)&(phir_min<2*pi-t2))
   Fs_min=-t1;
else
   Fs_min=min(Fs1,Fs2);
end
Hm_min=((rr/rm)*(Fpm+Fs_min)/(Rpm+Rg)-bpm)/(mu0*mur_pm);

% Compute maximum rotor flux density considering only radial flux
% at edges.  By symmetry, only points of
% positive magnetization (plus spacer) need to be considered.
% The candidate points are the PM transition locations (on either side)
% and at the point of peak stator MMF
B1=abs((Fpm+Fs1)/(Rpm+Rg));
B2=abs((Fpm+Fs2)/(Rpm+Rg));
B3=abs((Fs1)/(Rsp+Rg));
B4=abs((Fs2)/(Rsp+Rg));
phir_max=phii;
if (phir_max>2*pi)
   phir_max=phir_max-2*pi;
end
if (phir_max<0)
   phir_max=phir_max+2*pi;
end
if ((phir_max>pi)&(phir_max<2*pi))
   if ((phir_max>pi+t2)&(phir_max<2*pi-t2))
      B5=abs((Fpm+t1)/(Rpm+Rg));
   else
      B5=abs(t1/(Rsp+Rg));
   end
else
   phir_max=phir_max+pi; 
   if ((phir_max>pi+t2)&(phir_max<2*pi-t2))
      B5=abs((Fpm-t1)/(Rpm+Rg));
   else
      B5=abs(-t1/(Rsp+Rg));
   end
end
Brot_max1=max([B1 B2 B3 B4 B5]);

% Compute the rotor flux density manitude based on sectors aligned
% with stator teeth
Brot_rad=Phit*ns/(2*pi*rr*l);
Brot_tan=Phib/(l*(rr-ri));
Brot_max2=max(max(max(abs(Brot_rad))),max(max(abs(Brot_tan))));

% Choose between the two estimates
Brot_max=max(Brot_max1,Brot_max2);

%Maximum rotor tooth flux density
Bt_max=max(max(abs(Bt)));

%Maximum backiron flux density
Bb_max=max(max(abs(Bb)));

%Maximum tooth base flux density
Btb_max=max(max(abs(Btb)));

%Tooth core loss (hysteresis and eddy current)
[pBtqr,Btqr,qr] = B_versus_qr(Bt,phism_tc,qrm,P);
Ptcl=apt*ns*l*loss_density(Btqr,pBtqr,wr,SS);

%Backiron core loss (hysteresis and eddy current)
[pBbqr,Bbqr,qr] = B_versus_qr(Bb,phism_tc,qrm,P);
Pbcl=pi*(ro^2-rb^2)*l*loss_density(Bbqr,pBbqr,wr,SS);

% Plot the tooth flux density
if nargin==28

   figure(fn); 
   
   subplot(2,1,2);
   plot(qr*180/pi,Bbqr,'bx');
   xlabel('\theta_r, Degrees');
   ylabel('Backiron Flux Density, T')
   set(gca,'XTick',[0 60 120 180 240 300 360]);
   %axis([0 360 -2 2]);
   grid on;
   
   subplot(2,1,1);
   plot(qr*180/pi,Btqr,'bx');
   xlabel('\theta_r, Degrees');
   ylabel('Tooth Flux Density, T');
   set(gca,'XTick',[0 60 120 180 240 300 360]);
   temp=axis;
   %axis([0 360 temp(3:4)]);
   grid on;
   
end

end  % smpmsm_field_analysis.m                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         