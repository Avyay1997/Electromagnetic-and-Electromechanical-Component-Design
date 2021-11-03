function [BE] = BE_core(rh,wcp,ws,wop,hp,hb,MP,pmub)
% BE  Sets up data structure associated with a bobbin E core.
%
% [BE] = BE_core(rh,wcp,ws,wop,hp,hb,MP,pmub)
% [BE] = BE_core(rh,wcp,ws,wop,hp,hb,MP)
% [BE] = BE_core(rh,wcp,ws,wop,hp,hb)
%
% Inputs:
%  rh     = radius of central hole (m)
%  wcp    = width of center post (m)
%  ws     = width of slot (m)
%  wop    = width of outer post (m)
%  hp     = height of post (m)
%  hb     = height of base (m)
%  MP     = structure of material parameters (optional)
%  pmub   = pointer to permeability function (optional)
%
% Outputs:
% BE        = structure of bobin E-core parameters
%  BE.rh    = radius of central hole (m)
%  BE.wip   = width of center post (m)
%  BE.ws    = width of slot (m)
%  BE.wop   = width of outer post (m)
%  BE.hp    = height of post (m)
%  BE.hb    = height of base (m)
%  BE.rcp   = outside radius of center post (m)
%  BE.rop   = inside radius of outers posts (m)
%  BE.re    = smallest radius to outside edge (m)
%  BE.mcs_cp= magnetic cross section of center post (m^2);
%  BE.mcs_op= magnetic cross section of one edge post (m^2);
%  BE.mcs_b = magentic cross section of base (m^2);
%  BE.v     = volume of core (m^3)
%  BE.MP    = structure of material parameters (if MP supplied)
%  BE.m     = mass of core (kb) (if MP supplied)
%  BE.pmub  = pointer to permeability function (if supplied)
%
% Internal:
%
% Written by:
% S.D. Sudhoff
% School of Electrical and Computer Engineering
% 1285 Electrical Engineering Building
% West Lafayette, IN 47907-1285
% Phone: 765-494-3246
% Fax: 765-494-0676
% E-mail: sudhoff@ecn.purdue.edu

% assign known parameters
BE.rh=rh;
BE.wcp=wcp;
BE.ws=ws;
BE.wop=wop;
BE.hp=hp;
BE.hb=hb;

% compute radii
BE.rcp=rh+wcp;
BE.rop=BE.rcp+ws;
BE.re=BE.rop+wop;

% compute magnetic cross sections
BE.mcs_cp=pi*(BE.rcp^2 - BE.rh^2);
BE.mcs_op=2.0*BE.rcp*BE.wop;
BE.mcs_b =BE.hb*2*BE.rcp;

% compute volume
BE.v = (BE.mcs_cp + 2*BE.mcs_op)*BE.hp+ ...
       (4*BE.rcp*BE.re-pi*BE.rh^2)*BE.hb;   
   
% assign optional parameters and calculations
if nargin>6
   BE.MP=MP; 
   BE.m=BE.v*BE.MP.row;
end
if nargin==8
   BE.pmub=pmub;
end

%  Copyright 2013 - Scott Sudhoff 
% 
%  The program is distributed under the terms of the GNU Lesser
%  General Public License(GNU LGPL). 
% 
%  This file is part of MEC Toolbox 3.2
% 
%  MEC Toolbox is free software: you can redistribute it and/or modify
%  it under the terms of the GNU Lesser General Public License as published
%  by the Free Software Foundation, either version 3 of the License, or 
%  (at your option) any later version.
% 
%  MEC Toolbox is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
%  See the GNU Lesser General Public License for more details.
% 
%  You should have received a copy of the GNU Lesser General Public
%  License along with MEC Toolbox.  If not, see <http://www.gnu.org/licenses/>.
