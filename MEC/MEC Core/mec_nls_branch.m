function [MEC] = mec_nls_branch(MEC,bn,A,l,mi,Fs,PHIs,ml,PHIe)
% mec_nls_branch  Adds a magnetically non-linear explicit mesh branch to a 
%                 magnetic equivalent circuit
%
% [MEC]     = mec_nls_branch(MEC,bn,A,l,pmub,Fs,PHIs,ml,phi)
%
% Inputs:
% MEC      = MEC data structure (see mec_init for documentation)
% bn       = branch number.  This should be unique
% A        = effective cross sectional area of the branch
% l        = length of the branch
% mi       = material index
% Fs       = branch MMF source (A)
% PHIs     = branch flux source (Wb)
% ml       = mesh list for the branch
% PHIe     = estimate for the flux into branch (wb)
%
% Note: R and PHIs are in parallel, and this combination
%       is in series with Fs.  The postive node of Fs is
%       toward positive node of branch.  PHIs is directed
%       from the postive node towards the negative node.
%
% Outputs:
% MEC      = Updated MEC data structure
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
MEC.mb.t(bn)=2;
MEC.mb.R(bn)=NaN;
MEC.mb.Fs(bn)=Fs;
MEC.mb.PHIs(bn)=PHIs;
MEC.mb.A(bn)=A;
MEC.mb.l(bn)=l;
MEC.mb.mi(bn)=mi;
MEC.mb.PHIe(bn)=PHIe;

% find positive and negative branch lists
MEC.mb.pml{bn}=ml(ml>0);
MEC.mb.nml{bn}=-ml(ml<0);

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