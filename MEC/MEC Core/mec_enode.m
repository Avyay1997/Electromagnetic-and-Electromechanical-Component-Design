function [MEC] = mec_enode(MEC,n,ml)

% mec_enode  Adds an explicit node to a magnetic equivalent circuit
%
% [MEC]     = mec_enode(MEC,n,ml)
%
% Inputs:
% MEC      = MEC data structure (see mec_solve for documentation)
% n        = Node to add
% ml       = Mesh list.  This is a list of all explicit mesh fluxes
%            terminating or eminating from the node.  Explicit mesh
%            fluxes which terminate on the explicit node are considered 
%            positive.  Explicit mesh fluxes which eminate on the explicit
%            node are considered negative.  For example if mesh flux 5
%            terminates on a node, and mesh flux 7 eminates from a node,
%            then the mesh list is {5,-7}
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
MEC.nl.pml{n}= ml(ml>0);
MEC.nl.nml{n}=-ml(ml<0);

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