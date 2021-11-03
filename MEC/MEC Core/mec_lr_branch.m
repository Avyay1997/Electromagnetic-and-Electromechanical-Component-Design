function [MEC] = mec_lr_branch(MEC,bn,R,ml)
% mec_lr_branch  Adds a linear reluctance branch to a explicit mesh
%                list of a magnetic equivalent circuit
%
% [MEC]     = mec_lr_branch(MEC,bn,R,ml)
%
% Inputs:
% MEC      = MEC data structure (see mec_init for documentation)
% bn       = branch number.  This should be unique
% R        = branch reluctance (1/H)
% ml       = mesh list
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
MEC.mb.t(bn)=1;
MEC.mb.R(bn) =R;
MEC.mb.Fs(bn)=0;
MEC.mb.PHIs(bn)=0;

% find positive and negative branch lists
MEC.mb.pml{bn}= ml(ml>0);
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