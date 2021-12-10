function [G] = pmac_geometry(G)
% pmac_geometry performs a variety of geometry calculations for a permanent
%               magnet ac machine. Used in Chapter 9 of Power Magnetic 
%               Devices: A Multi-Objective Design Approach,"by S.D. Sudhoff
%
% Call:
% G=pmac_geometry(G)
%
% Inputs:
% M        = structure of perm. magnet parameters (see pmac_magnet_catalog)
% G        = a structure describing the machine geometry. input fields:
%  G.rrs   = shaft radius (m)
%  G.di    = depth of inert rotor (m)
%  G.drb   = depth of rotor backiron (m)
%  G.dtb   = depth of tooth base (m)
%  G.dm    = depth of magnet (m)
%  G.g     = air gap (m)
%  G.dttc  = depth of tooth tip at center (m)
%  G.dtte  = depth of tooth tip at edge (m)
%  G.alphat= tooth fraction (n/a)
%  G.alphatt tooth tip fraction (n/a)
%  G.dsb   = depth of stator backiron (m)
%  G.alphapm permanent magnet fraction (m)
%  G.l     = length of active stator structure (m)
%  G.P     = number of poles (integer)
%  G.Ss    = number of slots (integer)
%  G.phiss1= mechanical angle of position of center of slot 1 (rad)
%
% Outputs:
% G        = a structure describing the machine geometry. output fields:
%  G.phiss = slot center locations (mechanical) (rad)
%  G.phist = tooth center locations (mechanical) (rad)
%  G.rri   = radius rotor inert region (m)
%  G.rrb   = radius rotor backiron (m)
%  G.rrg   = radius rotor gap (m)
%  G.rst   = radius stator teeth (m)
%  G.rsi   = radius to stator inside tooth tip (m) 
%  G.qt    = angle of stator teeth (rad)
%  G.qtt   = angle of stator tooth tip (rad)
%  G.qst   = angle between tooth tips at tooth tip radius (rad)                              radius (rad)
%  G.wtb   = width of tooth base (m)
%  G.wtt   = width of tooth tip (m)
%  G.rsb   = radius of stator backiron (m)
%  G.qtb   = angle of tooth base (rad)
%  G.qti   = angle of tooth at inside tooth tip (rad) 
%  G.dst   = depth of stator tooth (m)
%  G.wso   = width of slot opening (m)
%  G.rss   = radius of stator shell (m)
%  G.aslt  = area of slot (m^2)
%  G.att   = area of tooth tip (m^2)
%  G.atb   = area of tooth base (m^2)
%  G.vst   = volume of stator teeth (m^3)
%  G.vsb   = volume of stator backiron (m^3)
%  G.vsl   = volume of stator laminatins (m^3)
%  G.vrb   = volume of rotor backiron (m^3)
%  G.rri   = volume of rotor inert region (m^3)
%  G.vpm   = volume of permanent magnet (m^3)
%  G.wttR  = width of tooth tip (RSA*) (m)
%  G.dttR  = depth of tooth tip (RSA*) (m)
%  G.wstR  = width of slot at tooth tip (RSA*) (m)
%  G.wsiR  = width of slot interior (RSA*) (m)
%  G.dsiR  = depth of slot interior (RSA*) (m)
%  G.wtbR  = width of tooth base (RSA*) (m)
%
%  * RSA = Rectangular slot approximation
%
% Internal:
% wct      = chord length across the top of slot
%            just under the tooth tip (m)
% wcb      = chord length across the bottom of slot(m)
% i        = Slot/tooth index
%
% Written by:
% S.D. Sudhoff
% School of Electrical and Computer Engineering
% 1285 Electrical Engineering Building
% West Lafayette, IN 47907-1285
% Phone: 765-494-3246
% Fax: 765-494-0676
% E-mail: sudhoff@ecn.purdue.edu

i=1:G.Ss;    
G.phiss =pi*(2*i-2)/G.Ss+G.phiss1;
G.phist =pi*(2*i-3)/G.Ss+G.phiss1;
G.rri=G.rrs+G.di;
G.rrb=G.rri+G.drb;
G.rrg=G.rrb+G.dm;
G.rst=G.rrg+G.g;
G.rsi=G.rst+G.dttc;
G.qt =2*pi*G.alphat/G.Ss;
G.qtt=2*pi*G.alphatt/G.Ss;
G.qst=2*pi/G.Ss-G.qtt;
G.wtb=2*G.rst*sin(G.qt/2);
G.wtt=2*G.rst*sin(G.qtt/2);
G.rsb=sqrt((G.wtb/2)^2+(G.rst*cos(G.qt/2)+G.dtb+G.dttc)^2);
G.qtb=2*asin(G.wtb/(2*G.rsb));   
G.qti=2*asin(G.wtb/(2*G.rsi)); 
G.dst=G.rsb-G.rst;
G.rss=G.rsb+G.dsb;
G.wso=2*G.rst*sin(G.qst/2);     
G.aslt=pi*(G.rsb^2-G.rsi^2)/G.Ss-G.wtb*G.dtb ...
       - G.rsb*(G.rsb*G.qtb-G.wtb*cos(G.qtb/2))/2 ...
       + G.rsi*(G.rsi*G.qti-G.wtb*cos(G.qti/2))/2;     
G.att=0.5*(2*G.wtt*G.dtte-G.rsi^2*G.qtt+G.wtt*G.rst*cos(G.qtt/2)+ ...
          (G.wtb+G.wtt)*(G.dttc-G.dtte)+ ...
          G.rsi^2*G.qti-G.wtb*G.rsi*cos(G.qti/2));
G.atb=pi*(G.rsb^2-G.rsi^2)/G.Ss-G.aslt;
G.vst=G.Ss*(G.att+G.atb)*G.l;
G.vsb=pi*(G.rss^2-G.rsb^2)*G.l;
G.vsl=G.vst+G.vsb;
G.vrb=pi*(G.rrb^2-G.rri^2)*G.l;                         
G.vri=pi*(G.rri^2-G.rrs^2)*G.l;                        
G.vpm=pi*(G.rrg^2-G.rrb^2)*G.alphapm*G.l;                
G.wttR=G.rst*G.qtt;
G.dttR=G.att/G.wttR;
wct=2.0*G.rsi*sin(pi/G.Ss-G.qti/2);   
wcb=2.0*G.rsb*sin(pi/G.Ss-G.qtb/2);   
G.wsiR=0.5*(wct+wcb);
G.wstR=G.rst*G.qst;
G.dsiR=G.aslt/G.wsiR;
G.wtbR=G.atb/G.dsiR;
  
end