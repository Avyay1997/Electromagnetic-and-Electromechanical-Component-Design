function [MEC] = mec_nl_material(MEC,mi,pmuf,mup)
% mec_nl_material  Adds a non-linear material to a magnetic equivalent
%                  circuit
%
% [MEC]     = mec_nl_material(MEC,mi,muf,mup)
%
% Inputs:
% MEC      = MEC data structure (see mec_init for documentation)
% mi       = material index.  Should be an integer from 1 to the number
%            of nonlinear materials in the circuit
% pmuf     = pointer to a function which calculated permeability and
%            the derivative of the permeability with respect to B as a
%            function of B
% mup      = parameters of permeability function
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

% test index
if (mi<1)||(mi>MEC.Nnlm)
   error('Invalid material index');
end
MEC.mat.muf{mi}=pmuf;
MEC.mat.mup{mi}=mup;


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
