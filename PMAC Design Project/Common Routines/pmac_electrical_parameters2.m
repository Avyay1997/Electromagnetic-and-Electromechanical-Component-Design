function [E] = pmac_electrical_parameters2(M,C,G,W)
% pmac_electrical_parameters computes the electrical parameters for a 
%                            surface mounted permanent magnet machine.
%                            Used in Chapter 9 of Power Magnetic Devices: A 
%                            Multi-Objective Design Approach," 
%                            by S.D. Sudhoff
%
% Call:
% E=pmac_electrical_parameters(M,C,G,W) 
%
% Inputs:
% M        = structure of perm. magnet parameters (see pmac_magnet_catalog)
% C        = structure of cond. information (see pmac_conductor_catalog)
% G        = structure of machine geometry (see pmac_geometry)
% W        = structure of winding variables (see pmac_winding)
%
% Outputs:
% E        = structure of lumped electrical parameters
%  E.Rs    = stator phase resistance (Ohms)        
%  E.Lq    = q-axis inductance (H)
%  E.Ld    = d-axis inductance (H)
%  E.lm    = d-axis flux linkage due to PM (Vs)
%
% Notes:
% On 12/2/2012 this version was modifed.  Formally, the
% smpmsm_lpm routine was called with G.dst as the argument
% for slot depth.  This is inconsisent with the use of the
% rectangular slot approximation.  This was thus replaced
% with G.dsiR
%
% Written by:
% S.D. Sudhoff
% School of Electrical and Computer Engineering
% 1285 Electrical Engineering Building
% West Lafayette, IN 47907-1285
% Phone: 765-494-3246
% Fax: 765-494-0676
% E-mail: sudhoff@ecn.purdue.edu

% Compute the lumped parameter model
if sum(abs(W.Nas))==0
   E.Rs=0.0;
   E.Lq=0.0;
   E.Ld=0.0;
   E.lm=0.0;
else
   [~,~,~,~,~,~,E.Lq,E.Ld,E.lm]= smpmsm_lpm( ...
                                       G.P,G.rrb,G.g,G.dm,G.dsb, ...
                                       G.dsiR,G.l,G.wtb,G.dttR, ...
                                       G.wttR,G.wstR,G.wsiR,W.dwR, ...
                                       C.sigma0,W.ac,W.leo,W.lew, ...
                                       G.alphapm,M.br,(1.0+M.x),1.0, ...                                            
                                       W.Ns1,W.Nas,W.Mas, ...
                                       W.Nbs,W.Mbs,G.rri);
   E.Rs=W.vcd/(C.sigma0*W.ac^2);
end