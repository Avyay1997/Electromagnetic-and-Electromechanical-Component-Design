function [MEC] = mec_lnb(MEC,n,pn,nn,P,Fs,PHIs)
% mec_lnb  Adds a magnetically-linear nodal branch to a 
%          magnetically equivalent circuit
%
% [MEC]     = mec_lmb(MEC,bn,pn,nn,P,Fs,PHIs)
%
% Inputs:
% MEC      = MEC data structure (see mec_solve for documentation)
% n        = branch number.  This should be unique for nodal branches
% pn       = positive explicit node for branch
% nn       = negative explicit node for branch
% P        = branch permeance (H)
% Fs       = branch MMF source (A)
% PHIs     = branch flux source (Wb)
%
% Note: P and PHIs are in parallel, and this combination
%       is in series with Fs.  The postive node of Fs is
%       toward positive node of branch.  PHIs is directed
%       from the postive node towards the negative node.
%
% Outputs:
% MEC      = Updated MEC data structure (see mec_solve for documentation)
%
% Written by:
% S.D. Sudhoff
% School of Electrical and Computer Engineering
% 1285 Electrical Engineering Building
% West Lafayette, IN 47907-1285
% Phone: 765-494-3246
% Fax: 765-494-0676
% E-mail: sudhoff@ecn.purdue.edu

% update appropriate branch fields
MEC.nb.P(n) =P;
MEC.nb.Fs(n)=Fs;
MEC.nb.PHIs(n)=PHIs;
MEC.nb.pn(n)=pn;
MEC.nb.nn(n)=nn;

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
%  License along with MEC Toolbox.  
%  If not, see <http://www.gnu.org/licenses/>.