function [F] = pmac_field_analysis2(M,G,W,IIs,Iphii,S,J,wrm)
% pmac_field_analysis analyzes the fields in a surface mounted permanent
%                     magnet ac machine. Used in Chapter 9 of  "Power 
%                     Magnetic Devices: A Multi-Objective Design Approach,"
%                     by S.D. Sudhoff
%
% Call:
% F=pmac_field_analysis(M,G,W,I,S,J,wrm)
% 
% Inputs:
% M        = structure of perm. magnet parameters (see pmac_magnet_catalog)
% G        = structure of machine geometry (see pmac_geometry)
% W        = structure of winding variables (see pmac_winding)
% I        = structure of current informatio (see pmac_current)
% S        = structure of steel parameters (see pmac_steel_catalog)
% J        = number of rotor positions to use
% wrm      = mechanical rotor speed (rad/s)
%
% Outputs:
% F        = structure of field information
%  F.qrc   = electrical rotor position over a cycle (rad)
%  F.Bt1c  = flux density in tooth 1 over a cycle (T)
%  F.pBt1c = derivative of flux density in tooth 1 with
%            respect to electrical rotor position (T/s)
%  F.Bb1c  = flux density in backiron segment 1 over a cycle (T)
%  F.pBb1c = derivative of flux density in tooth 1 with
%            respect to electrical rotor position (T/s)
%  F.Brbtmx= maximum tangential flux density in rotor backiron (T)
%  F.Brbrmx= maximum radial flux density in rotor backrion (T)
%  F.Hmn   = minimum field intensity in positively magnetized region of 
%            magnet (A/m)
%  F.Pc    = core loss (W)
%
% Internal:
% ind      = index of first half cycle
% wr       = electrical rotor position (rad/s)
% pt       = core loss density in teeth (W/m^3)
% pb       = core loss density in backiron (W/m^3)
%
% Written by:
% S.D. Sudhoff
% School of Electrical and Computer Engineering
% 1285 Electrical Engineering Building
% West Lafayette, IN 47907-1285
% Phone: 765-494-3246
% Fax: 765-494-0676
% E-mail: sudhoff@ecn.purdue.edu

[F.qrc,F.Bt1c,F.pBt1c,F.Bb1c,F.pBb1c, ...
          F.Brbtmx,F.Brbrmx,F.Hmn]= ...
smpmsm_field_analysis4(J, ...
                       G.P,G.Ss,G.l,G.rrb,G.rri,G.rrg, ...
                       G.wtb,G.dsb, ...
                       0.0,W.Ns1,G.dm,G.alphapm,M.br, ...
                       (1+M.x),1.0,G.g, ...
                       IIs,Iphii, ...
                       G.phist);
if (wrm==0)
   F.Pc=0;
else
   ind=[1:length(F.Bt1c)/2+1]; 
   wr=wrm*G.P/2;
   pt=core_loss_density(F.qrc(ind),F.Bt1c(ind),F.pBt1c(ind)*wr,wr,S);
   pb=core_loss_density(F.qrc(ind),F.Bb1c(ind),F.pBb1c(ind)*wr,wr,S);
   F.Pc=pt*G.vst+pb*G.vsb;
%    disp('old')
%    [pt pb]
%    figure(102)
%    plot(F.qrc(ind),F.Bt1c(ind))
%    figure(103)
%    plot(F.qrc(ind),F.pBt1c(ind)*wr)
end

end

function pld=core_loss_density(qe,B,dBdt,we,CM)
% core_loss_density    computes the core loss density of a material
%                      based on a 1/2 cycle data. Assumes B is symmetric
%                      and half-wave symmetric.
%
% pld=core_loss_density(qe,B,dBdt,we,CM)
%
% Inputs:
% qe         = angle of data ranging over 1/2 of a cycle (rad)
%              Should cover a full half cycle (last point = first point
%              plus 1/2 cycle)
% B          = flux density structure with data over 1/2 of a cycle (T)
% dBdt       = time derivative of flux density waveform (T/s)
% we         = radian frequency of data
% CM         = stucture of material data (see pmac_steel_catalog)
% 
% Outputs:
% pld        = power loss density (W/m^3)
% 
% Internal:
% f          = fundamental frequency (Hz)
% t          = time (s) (spans 1/4 cycle)
% N          = number of points in vectors
% i1,i2      = index arrays
% pB         = time derivative of B (T/s)
% pB2        = square of time derivative of B (T/s)^2
% Bmx        = maximum flux density (T)
% Bpk        = peak flux density (T)
% DB         = peak-to-peak flux density range (T)
% int_dBdt2  = integral of square of derivative of B over a cycle
% feq        = equivalent frequency
%
% Written by:
% S.D. Sudhoff
% School of Electrical and Computer Engineering
% 1285 Electrical Engineering Building
% West Lafayette, IN 47907-1285
% Phone: 765-494-3246
% Fax: 765-494-0676
% E-mail: sudhoff@ecn.purdue.edu

% Determine fundamental frequency
f=abs(we)/(2.0*pi);
   
% Determine elapsed time
t=abs(qe/we);
t=t-t(1);
 
% Determine the square of the time derivative of flux density
pB2=dBdt.^2;                     
   
% Determine maximum B, minimum B, peak B, and range of B over a cycle
Bmx=max(abs(B));
Bpk=Bmx;
DB=2.0*Bmx;
    
% Determine integral (factor of 2 is because using 1/2 of a cycle)
int_dBdt2=2.0*trapz(t,pB2);
  
% Determine equivalent frequency
feq=2*int_dBdt2/(DB*pi).^2;
   
% Determine output
pld=CM.MSE.kh*Bpk^(CM.MSE.beta)*feq^(CM.MSE.alpha-1)*f+ ...
    CM.MSE.ke*f*int_dBdt2;

end


