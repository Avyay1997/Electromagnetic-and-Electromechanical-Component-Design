function Wf = mec_field_energy(PHImb,Fmb,PHInb,Fnb,MEC)
% mec_field_energy calculates field energy in a MEC
%                  Assumptions:
%                  1.) The MMF source in the nodal and mesh branches
%                      represents an external current and is not used
%                      to represent magnetic material
%                  2.) The Flux source in the mesh branches represents
%                      permanent magnet material
%                  3.) The Flux source in the nodal branches represents
%                      permanent magnet material UNLESS the permeance
%                      for the branch is zero.
%
% Inputs:
% MEC      = MEC data structure (see mec_init for description of fields)
% PHImb    = Vector of explicit mesh branch fluxes (MEC.Nemb by 1) (Wb)
% Fmb      = Vector of explicit mesh branch MMF drops (MEC.Nemb by 1) (At)
% Fnb      = Vector of node branch MMF drops (At)
% PHInb    = Vector of fluxes into nodal branches (Wb)
% 
% Outputs:
% Wf       = Field energy (J)
% 
% Intenal:
% Em       = Energy stored in mesh branches
% index    = Nodal branches with non zero permeance
% En       = Energy stored in nodal branches

   Em=0.5*PHImb'*(Fmb-MEC.mb.Fs);
   index=MEC.nb.P>0;
   En=0.5*PHInb(index)'*(Fnb(index)-MEC.nb.Fs(index));
   Wf=Em+En;
%  keyboard

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
%  License along with MEC Toolbox.  
%  If not, see <http://www.gnu.org/licenses/>.