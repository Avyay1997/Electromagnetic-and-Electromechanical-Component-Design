function [ts] = mec_test(MEC)
% mec_test   Performs some tests on the MEC setup.
%
% ts = mec_test(MEC)
%
% Inputs:
% MEC      = Magnetic equivalent circuit structure (see mec_init for doc)
%
% Output:
% ts       = Test status.  One if test passed; zero if not.
%
% Internal:
% bn       = Branch number
% nn       = Node number
% mn       = Mesh number
% ml       = Mesh list
% tf1...   = Test flag.  High (1) if test failed.
% tf20
%
% Written by:
% S.D. Sudhoff
% School of Electrical and Computer Engineering
% 1285 Electrical Engineering Building
% West Lafayette, IN 47907-1285
% Phone: 765-494-3246
% Fax: 765-494-0676
% E-mail: sudhoff@ecn.purdue.edu

% Initialize test status to passed
ts=1;

% Test mesh branches
for bn=1:MEC.Nemb
   tf1 =isnan(MEC.mb.R(bn));
   tf2 =isnan(MEC.mb.Fs(bn));
   tf3 =isnan(MEC.mb.PHIs(bn));
   tf4 =isnan(MEC.mb.A(bn));
   tf5 =isnan(MEC.mb.l(bn));
   tf6 =(MEC.mb.mi<1)|(MEC.mb.mi>MEC.Nnlm);
   tf7 =0;
   tf8 =isnan(MEC.mb.PHIe(bn));
   ml=[MEC.mb.pml{bn} MEC.mb.nml{bn}];
   tf9 =isempty(ml);
   tf10 = (MEC.mb.mi(bn)>MEC.Nnlm)|(MEC.mb.mi(bn)<0);
   if ~tf9
      tf9=~((min(ml)>=1)&(max(ml)<=MEC.Nem));
   end
   if tf2
      disp(['Mesh branch ' num2str(bn) ' MMF source has not been specified']);             
   end
   if tf3
      disp(['Mesh branch ' num2str(bn) ' flux source has not been specified']);             
   end 
   switch MEC.mb.t(bn)
      case 1,
         if tf1 
            disp(['Mesh branch ' num2str(bn) ' reluctance has not been specified']);
         end
         if ~(tf5&&tf6&&tf7&&tf8)
            disp(['Mesh branch ' num2str(bn) ' is linear branch with nonlinear fields specified']); 
         end
         ts=ts&(~tf1)&(~tf2)&(~tf3)&(tf5&&tf6&&tf7&&tf8);
      case 2,
         if tf4
            ts=0;
            disp(['Mesh branch ' num2str(bn) ' magnetic cross section has not been specified']);
         end   
         if tf5
            ts=0;
            disp(['Mesh branch ' num2str(bn) ' magnetic length has not been specified']);
         end  
         if tf6
            ts=0;
            disp(['Mesh branch ' num2str(bn) ' permeability function has not been specified']);
         end  
         if tf7
            ts=0;
            disp(['Mesh branch ' num2str(bn) ' permeability parameters have not been specified']);
         end
         if tf8
            ts=0;
            disp(['Mesh branch ' num2str(bn) ' branch flux estimate has not been specified']);
         end  
         if ~(tf1)
            ts=0;
            disp(['Mesh branch ' num2str(bn) ' is nonlinear branch with linear field specified']); 
         end
         ts=ts&(~tf4)&(~tf5)&(~tf6)&(~tf7)&(~tf8)&(tf1);
      otherwise
         ts=0;
         disp(['Mesh branch ' num2str(bn) ' has invalid type. ']);
   end
   if tf9
      disp(['Mesh branch ' num2str(bn) ' has an invalid mesh list.']);
   end
   if tf10
      disp(['Mesh branch ' num2str(bn) ' has invalid nonlinear material.']);
   end
   ts=ts&(~tf9)&(~tf10);

end

% Test nodal branches
for bn=1:MEC.Nenb
   tf11=isnan(MEC.nb.pn(bn));
   if (~tf11)
      tf11=~((MEC.nb.pn(bn)>=0)&(MEC.nb.pn(bn)<=MEC.Nen));
   end
   tf12=isnan(MEC.nb.nn(bn));
   if (~tf12)
      tf12=~((MEC.nb.nn(bn)>=0)&(MEC.nb.nn(bn)<=MEC.Nen));
   end  
   tf13=isnan(MEC.nb.P(bn));
   tf14=isnan(MEC.nb.Fs(bn));
   tf15=isnan(MEC.nb.PHIs(bn));
   if tf11
      disp(['Nodal branch ' num2str(bn) ' positive node is not valid.']);
   end
   if tf12
      disp(['Nodal branch ' num2str(bn) ' negative node is not valid.']);
   end   
   if tf13
      disp(['Nodal branch ' num2str(bn) ' permeance not specified.']);
   end 
   if tf14
      disp(['Nodal branch ' num2str(bn) ' MMF source not specified.']);
   end 
   if tf15
      disp(['Nodal branch ' num2str(bn) ' Flux source not specified.']);
   end 
   ts=ts&(~tf10)&(~tf11)&(~tf12)&(~tf13)&(~tf14)&(~tf15);
end

% Test mesh list to make sure all mesh fluxes appear
% in at least one mesh branch
for mn=1:MEC.Nem
   tf16=1;
   for bn=1:MEC.Nemb
      if ~isempty(intersect(mn,MEC.mb.pml{bn}))
         tf16=0;
      end
      if ~isempty(intersect(mn,MEC.mb.nml{bn}))
         tf16=0;
      end
   end
   if tf16
      disp(['Mesh flux ' num2str(mn) ' does not appear in any mesh branch.']);
   end
   ts=ts&(~tf16);   
end

% Test node list to make sure all nodes appear in at least
% one nodal branch
for nn=1:MEC.Nen
   ml=[MEC.nl.pml{nn} MEC.nl.nml{nn}];
   tf17=isempty(ml);
   if ~tf17
      tf17=~((min(ml)>=1)&(max(ml)<=MEC.Nem));
   end
   if tf17
      disp(['Node ' num2str(nn) ' has invalid mesh list.']);
   end 
   ts=ts&(~tf17);
   tf18=1;
   for bn=1:MEC.Nenb
      if MEC.nb.pn(bn)==nn;
         tf18=0;
      end
      if MEC.nb.nn(bn)==nn;
         tf18=0;
      end
   end
   if tf18
      disp(['Node ' num2str(nn) ' is not used in any nodal branch.']);
   end    
   ts=ts&(~tf18);
end

% Test material list
for mi=1:MEC.Nnlm
   tf19=isempty(MEC.mat.muf{mi});
   if tf19
      disp(['Permeabilty function for material ' num2str(mi) ...
            ' not specified']);
   end
   ts=ts&(~tf19);
   tf20=isempty(MEC.mat.mup{mi});
   if tf20
      disp(['Permeability parameters for material ' num2str(mi) ...
            ' not specified']);
   end
   ts=ts&(~tf20);
end
   
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