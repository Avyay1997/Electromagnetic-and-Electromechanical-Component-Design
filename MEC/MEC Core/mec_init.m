function [MEC] = mec_init(varargin)
% mec_init   Initializes an Magnetic Equivalent Circuit analysis
%            of a magnetic structure.  All fields are initialized to zero
%
% [MEC]     = mec_init(Nem,Nemb)
% [MEC]     = mec_init(Nem,Nemb,Nnlm)
% [MEC]     = mec_init(Nem,Nemb,Nen,Nenb,Nnlm)
% [MEC]     = mec_init(Nem,Nemb,Nen,Nenb,Npsb,Nnlm)
%
% Inputs:
% Nem      = Number of explicit meshes
% Nemb     = Number of explicit mesh branches
% Nen      = Number of explicit nodes
% Nenb     = Number of explicit node branches
% Npsb     = Number of pole symmetric branches
% Nnlm     = Number of different non-linear magentic materials
%
% Output:
% MEC      = Initialized MEC data structure
%   MEC.Nem     = Number of explicit meshes
%   MEC.Nemb    = Number of explicit mesh branches
%   MEC.Nen     = Number of explicit nodes
%   MEC.Nenb    = Number of explicit node branches
%   MEC.mxit    = Maximum number of iterations
%   MEC.rec     = Relative error criteria
%   MEC.aec     = Absolute error criteria
%   MEC.mat.*     Non-linear magnetic materials structure
%   MEC.mat.muf = Cell array with pointer to permeability functions
%                 (mu and the partial of mu with respect to B as
%                 functions of B)
%   MEC.mat.mup = Cell array of structures with parameters for the
%                 permeabiity function
%   MEC.mb.*      Mesh branch structure.  Made of vectors associated
%                 with the explicit mesh branch list (members Nemb by 1)
%   MEC.mb.t    = Mesh branch types.  1=linear; 2=non-linear
%   MEC.mb.R    = Mesh branch reluctances (1/H)
%   MEC.mb.Fs   = Mesh branch MMF sources (At)
%   MEC.mb.PHIs = Mesh branch flux sources (Wb)
%   MEC.mb.A    = Mesh branch magnetic cross sections (m^2)
%   MEC.mb.l    = Mesh branch magnetic lengths (m)
%   MEC.mb.pmub = Mesh branch pointers to permeability functions 
%   MEC.mb.MMP  = Mesh branch magnetic material parameters 
%   MEC.mb.PHIe = Estimated flux into branches (Wb)
%   MEC.mb.pml  = Cell vector with postive mesh list for branch
%   MEC.mb.nml  = Cell vector with negative mesh list for branch
%   MEC.mb.mi   = material index for branch (irrelevant if linear material)
%   MEC.nb.*      Node branch structure. Made up of vectors associated 
%                 with the explicit node branch list (members Nenb by 1)
%   MEC.nb.pn   = Node branch positive node
%   MEC.nb.nn   = Node branch negative node
%   MEC.nb.P    = Node branch permeance (H)
%   MEC.nb.Fs   = Node branch MMF source (A)
%   MEC.nb.PHIs = Node branch flux source (Wb)
%   MEC.psb.*     Pole symmetric branch structure. Made up of vectors  
%                 associated with the explicit node pole symmetric list 
%                 (members Npsb by 1)
%   MEC.psb.pn1 = Pole symmetric branch positive node 1
%   MEC.psb.pn2 = Pole symmetric branch positive node 2
%   MEC.psb.P   = Pole symmetric branch permeance (H)
%   MEC.psb.Fs  = Pole symmetric branch MMF source (A)
%   MEC.psb.PHIs= Pole symmetric branch flux source (Wb)
%   MEC.nl.*      Explicit node list strucure. Members are (Nen by 1)
%   MEC.nl.pml  = Positive mesh list
%   MEC.nl.nml  = Negative mesh list
%
% Written by:
% S.D. Sudhoff
% School of Electrical and Computer Engineering
% 1285 Electrical Engineering Building
% West Lafayette, IN 47907-1285
% Phone: 765-494-3246
% Fax: 765-494-0676
% E-mail: sudhoff@ecn.purdue.edu

% Solver parameters
switch nargin
   case 2,
      MEC.Nem  = varargin{1};
      MEC.Nemb = varargin{2};
      MEC.Nen  = 0;
      MEC.Nenb = 0;
      MEC.Npsb = 0;
      MEC.Nnlm = 0;
   case 3,
      MEC.Nem  = varargin{1};
      MEC.Nemb = varargin{2};
      MEC.Nen  = 0;
      MEC.Nenb = 0;
      MEC.Npsb = 0;
      MEC.Nnlm = varargin{3};
   case 5,
      MEC.Nem  = varargin{1};
      MEC.Nemb = varargin{2};
      MEC.Nen  = varargin{3};
      MEC.Nenb = varargin{4};
      MEC.Npsb = 0;
      MEC.Nnlm = varargin{5};       
   case 6,
      MEC.Nem  = varargin{1};
      MEC.Nemb = varargin{2};
      MEC.Nen  = varargin{3};
      MEC.Nenb = varargin{4};
      MEC.Npsb = varargin{5};
      MEC.Nnlm = varargin{6};
   otherwise
      error('Inappropriate number of input arguments');
end % switch

% Solver tolerances
MEC.mxit=50;
MEC.rec=1e-12;
MEC.aec=1e-14;

% Initialize fields associated with non-linear materials list
MEC.mat.muf = cell(MEC.Nnlm,1);
MEC.mat.mup = cell(MEC.Nnlm,1);

% Initialize fields associated with explicit mesh branches
MEC.mb.t   = NaN(MEC.Nemb,1);
MEC.mb.R   = NaN(MEC.Nemb,1);
MEC.mb.Fs  = NaN(MEC.Nemb,1);
MEC.mb.PHIs= NaN(MEC.Nemb,1);
MEC.mb.A   = NaN(MEC.Nemb,1);
MEC.mb.l   = NaN(MEC.Nemb,1);
MEC.mb.pmub= cell(MEC.Nemb,1);
MEC.mb.MMP = cell(MEC.Nemb,1);
MEC.mb.PHIe= NaN(MEC.Nemb,1);
MEC.mb.pml = cell(MEC.Nemb,1);
MEC.mb.nml = cell(MEC.Nemb,1);
MEC.mb.mi  = NaN(MEC.Nemb,1);

% Initialize fields associated with explicit node branches
MEC.nb.pn  = NaN(MEC.Nenb,1);
MEC.nb.nn  = NaN(MEC.Nenb,1);
MEC.nb.P   = NaN(MEC.Nenb,1);
MEC.nb.Fs  = NaN(MEC.Nenb,1);
MEC.nb.PHIs= NaN(MEC.Nenb,1);

% Initialize fields associated with pole symmetric branches
MEC.psb.pn1 = NaN(MEC.Npsb,1);
MEC.psb.pn2 = NaN(MEC.Npsb,1);
MEC.psb.P   = NaN(MEC.Npsb,1);
MEC.psb.Fs  = NaN(MEC.Npsb,1);
MEC.psb.PHIs= NaN(MEC.Npsb,1);

% Initialize fields assoicated with explicit node lists
MEC.nl.pml = cell(MEC.Nen,1);
MEC.nl.nml = cell(MEC.Nen,1);

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