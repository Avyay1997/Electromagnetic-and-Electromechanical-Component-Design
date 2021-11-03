function [varargout] = mec_solve(MEC)
% mec_solve  Solves the magnetic equivalent circuit
%
% [PHIm,PHImb,PHImR,Fmb,c] = mec_solve(MEC)
% [PHIm,PHImb,PHImR,Fmb,Fn,Fnb,PHInb,PHInP,c] = mec_solve(MEC)
% [Fn,Fnb,PHInb,c] = mec_solve(MEC)
%
% Inputs:
% MEC      = MEC data structure (see mec_init for description of fields)
%
% Outputs:
% PHIm     = Vector of explicit mesh fluxes (MEC.Nem by 1) (Wb)
% PHImb    = Vector of explicit mesh branch fluxes (MEC.Nemb by 1) (Wb)
% PHImR    = Vector of explicit mesh reluctance fluxes (MEC.Nemb by 1)(Wb)
% Fmb      = Vector of explicit mesh branch MMF drops (MEC.Nemb by 1) (At)
% Fn       = Vector of explicit node potentials (MEC.Nen by 1) (At)
% Fnb      = Vector of node branch MMF drops (At)
% PHInb    = Vector of fluxes into nodal branches (Wb)
% PHInP    = Vector of fluxes into nodal branch permeances (Wb)
% c        = Convergence. 1 if analysis converges; 0 otherwise
%
% Internal:
% A        = System matrix (used to express MEC system equations)
% b        = System forcing vector (used to express MEC system equations)
% x        = Solution of system equations
% pml      = Positive mesh list for branch
% nml      = Negative mesh list for branch
% bn       = Branch number
% R        = Reluctance of a branch (H^-1)
% st       = Composite source term for a branch (A)
% pn       = Positive node
% pni      = Positive node index
% nn       = Negative node
% nni      = Negative node index
% n        = Explicit node
% ni       = Node index
% A0       = magnetically linear part of system matrix (H^-1)
% b0       = magnetically linear part of source vector (A)
% PHImht   = previous estimate of PHIm (Wb)
% i        = index variable
% k        = iteration count
% Ab       = cross sectional area of a magnetically non-linear branch (m^2)
% lb       = effective length of a magnetically non-linear branch (m)
% PHI0    =  flux through a magnetically non-linear branch 
%            at point of linearization (Wb)
% B0       = flux density in a magnetically non-linear branch at point
%            of linearization (T)
% mu       = permeability of magnetically non-linear branch at B0 (H/m)
% pmu      = derivative of permeability with respect to flux density 
%            at B0 (H/(Tm))
% Abmu     = Product of Ab amd mu (m H)
% R0       = Absolute (not incremental) reluctance at B0  (H^-1)
% S0       = Temporary term at B0
% Fsmbeff  = Vector of effective source MMF (branch is viwed as a 
%            simple Thevenin equivalent) (A t)
% Rmbeff   = Vector of incremental relutances of branches
% e        = a measure of the solution error
% emlb     = set of explicit mesh linear branches
% emnb     = set of explicit mesh non-linear branches
%
% Written by:
% S.D. Sudhoff
% School of Electrical and Computer Engineering
% 1285 Electrical Engineering Building
% West Lafayette, IN 47907-1285
% Phone: 765-494-3246
% Fax: 765-494-0676
% E-mail: sudhoff@ecn.purdue.edu

% Initialize matrices
A=zeros(MEC.Nem+MEC.Nen,MEC.Nem+MEC.Nen);
b=zeros(MEC.Nem+MEC.Nen,1);
if (MEC.Nemb>0)
   PHIm(MEC.Nem,1)=0;
end
if (MEC.Nemb>0)
   PHImb(MEC.Nemb,1)=0;
   PHImR(MEC.Nemb,1)=0;
   Fmb(MEC.Nemb,1)=0;
end
if (MEC.Nen>0)
   Fn(MEC.Nen,1)=0;
end
if (MEC.Nenb>0)
   Fnb(MEC.Nenb,1)=0;
end
Fsmbeff=MEC.mb.Fs-MEC.mb.R.*MEC.mb.PHIs;
Rmbeff=MEC.mb.R;
   
% Divide problem by type
if max(MEC.mb.t)==1

   % Magnetically Linear Problem-------------------------------------------
   
   % Branches associated with explicit mesh
   if (MEC.Nemb>0)
      for bn=1:MEC.Nemb
       
         % get the mesh lists
         pml=MEC.mb.pml{bn};
         nml=MEC.mb.nml{bn};
        
         % update reluctance matrix
         R=MEC.mb.R(bn);
         A(pml,pml)=A(pml,pml)+R;
         A(pml,nml)=A(pml,nml)-R;
         A(nml,pml)=A(nml,pml)-R;
         A(nml,nml)=A(nml,nml)+R;

         % update source matrix
         st=R*MEC.mb.PHIs(bn)-MEC.mb.Fs(bn);
         b(pml)=b(pml)+st;
         b(nml)=b(nml)-st;
    
      end   
   end
 
   % Branches associated with explicit nodes
   if (MEC.Nenb>0)
      for bn=1:MEC.Nenb
          
          % get the positive and negative nodes and matrix indices
          pn=MEC.nb.pn(bn);
          nn=MEC.nb.nn(bn);
          pni=pn+MEC.Nem;
          nni=nn+MEC.Nem;
          
          % update permeance and source matrix
          P=MEC.nb.P(bn);
          st=P*MEC.nb.Fs(bn)-MEC.nb.PHIs(bn);
          if (pn>0)
             A(pni,pni)=A(pni,pni)+P;
             b(pni)=b(pni)+st;
          end
          if (nn>0)
             A(nni,nni)=A(nni,nni)+P; 
             b(nni)=b(nni)-st;
          end
          if (pn>0)&&(nn>0)
             A(pni,nni)=A(pni,nni)-P;
             A(nni,pni)=A(nni,pni)-P;
          end
                       
      end
   end
   
   % Pole symmetric branches of explicit nodes
   if (MEC.Npsb>0)
      for bn=1:MEC.Npsb
          
          % get the nodes and matrix indices
          pn1=MEC.psb.pn1(bn);
          pn2=MEC.psb.pn2(bn);
          pn1i=pn1+MEC.Nem;
          pn2i=pn2+MEC.Nem;
          
          % update permeance and source matrix
          P=MEC.psb.P(bn);
          st=P*MEC.psb.Fs(bn)-MEC.psb.PHIs(bn);
          if (pn1>0)
             A(pn1i,pn1i)=A(pn1i,pn1i)+P;
             b(pn1i)=b(pn1i)+st;
          end
          if (pn2>0)
             A(pn2i,pn2i)=A(pn2i,pn2i)+P; 
             b(pn2i)=b(pn2i)+st;
          end
          if (pn1>0)&&(pn2>0)
             A(pn1i,pn2i)=A(pn1i,pn2i)+P;
             A(pn2i,pn1i)=A(pn2i,pn1i)+P;
          end
                       
      end
   end
   
   
   % Explicit-mesh, explicit node interaction
   if (MEC.Nenb>0)&&(MEC.Nemb>0)
      for n=1:MEC.Nen
          pml=MEC.nl.pml{n};
          nml=MEC.nl.nml{n};
          ni=MEC.Nem+n;
          A(pml,ni)=A(pml,ni)+1;
          A(nml,ni)=A(nml,ni)-1;
          A(ni,pml)=A(ni,pml)-1;
          A(ni,nml)=A(ni,nml)+1;
      end
   end
   
   % Solve system equations
   x=A\b;
   PHIm=x(1:MEC.Nem);
   Fn=x(MEC.Nem+1:end);

   % determine the branch fluxes
   for bn=1:MEC.Nemb
       pml=MEC.mb.pml{bn};
       nml=MEC.mb.nml{bn};
       PHImb(bn)=sum(PHIm(pml))-sum(PHIm(nml));
   end
      
   % linear case always coverges
   c=1;
   
   % END Magnetically Linear Problem---------------------------------------
        
else
    
   % Magnetically Non-Linear Problem---------------------------------------
   
   % Step 0 - Initialize extra matrice/vectors for nonlinear case
   A0=zeros(MEC.Nem+MEC.Nen,MEC.Nem+MEC.Nen);
   b0=zeros(MEC.Nem+MEC.Nen,1);
  
   % Step 1a - process magnetically linear mesh branches
   emlb=find(MEC.mb.t==1);
   for i=1:length(emlb)
       
       % determine branch numbers and lists
       bn=emlb(i);
       pml=MEC.mb.pml{bn};
       nml=MEC.mb.nml{bn};
        
       % update system matrix
       R=MEC.mb.R(bn);
       A0(pml,pml)=A0(pml,pml)+R;
       A0(pml,nml)=A0(pml,nml)-R;
       A0(nml,pml)=A0(nml,pml)-R;
       A0(nml,nml)=A0(nml,nml)+R;

       % update source vector
       st=R*MEC.mb.PHIs(bn)-MEC.mb.Fs(bn);
       b0(pml)=b0(pml)+st;
       b0(nml)=b0(nml)-st;
     
   end
   
   % Steb 1b - process explicit nodal branches
   if (MEC.Nenb>0)
      for bn=1:MEC.Nenb
          
          % get the positive and negative nodes and matrix indices
          pn=MEC.nb.pn(bn);
          nn=MEC.nb.nn(bn);
          pni=pn+MEC.Nem;
          nni=nn+MEC.Nem;
          
          % update permeance and source matrix
          P=MEC.nb.P(bn);
          st=P*MEC.nb.Fs(bn)-MEC.nb.PHIs(bn);
          if (pn>0)
             A0(pni,pni)=A0(pni,pni)+P;
             b0(pni)=b0(pni)+st;
          end
          if (nn>0)
             A0(nni,nni)=A0(nni,nni)+P; 
             b0(nni)=b0(nni)-st;
          end
          if (pn>0)&&(nn>0)
             A0(pni,nni)=A0(pni,nni)-P;
             A0(nni,pni)=A0(nni,pni)-P;
          end
                       
      end
   end

   % Step 1c
   % Pole symmetric branches of explicit nodes
   if (MEC.Npsb>0)
      for bn=1:MEC.Npsb
          
          % get the nodes and matrix indices
          pn1=MEC.psb.pn1(bn);
          pn2=MEC.psb.pn2(bn);
          pn1i=pn1+MEC.Nem;
          pn2i=pn2+MEC.Nem;
          
          % update permeance and source matrix
          P=MEC.psb.P(bn);
          st=P*MEC.psb.Fs(bn)-MEC.psb.PHIs(bn);
          if (pn1>0)
             A0(pn1i,pn1i)=A0(pn1i,pn1i)+P;
             b0(pn1i)=b0(pn1i)+st;
          end
          if (pn2>0)
             A0(pn2i,pn2i)=A0(pn2i,pn2i)+P; 
             b0(pn2i)=b0(pn2i)+st;
          end
          if (pn1>0)&&(pn2>0)
             A0(pn1i,pn2i)=A0(pn1i,pn2i)+P;
             A0(pn2i,pn1i)=A0(pn2i,pn1i)+P;
          end
                       
      end
   end
   
      
   % Steb 1d - explicit-mesh, explicit node interaction
   % This is also a linear term
   if (MEC.Nenb>0)&&(MEC.Nemb>0)
      for n=1:MEC.Nen
          pml=MEC.nl.pml{n};
          nml=MEC.nl.nml{n};
          ni=MEC.Nem+n;
          A0(pml,ni)=A0(pml,ni)+1;
          A0(nml,ni)=A0(nml,ni)-1;
          A0(ni,pml)=A0(ni,pml)-1;
          A0(ni,nml)=A0(ni,nml)+1;
      end
   end   
   
   % Step 2 - iteratively solve the system
   % initialize iteration count and error
   k=0;
   e=0;
   emnb=find(MEC.mb.t==2);
   PHImb(emnb)=MEC.mb.PHIe(emnb);
   while (k<=1)||((k<=MEC.mxit)&&(e>0))
       
      % Copy of initial solution
      PHImht=PHIm;
             
      % prepare system matrices
      A=A0;
      b=b0;
      
      % linearize the magnetics one material at a time
      for mi=1:MEC.Nnlm
          index=MEC.mb.mi==mi;
          Ab=MEC.mb.A(index);
          lb=MEC.mb.l(index);
          PHI0=PHImb(index);
          B0=PHI0./Ab;
          [mu,pmu]=feval(MEC.mat.muf{mi},MEC.mat.mup{mi},B0);
          Abmu=Ab.*mu;
          R0=lb./Abmu;
          S0=pmu.*(PHI0-MEC.mb.PHIs(index))./Abmu;
          Fsmbeff(index)=MEC.mb.Fs(index)+R0.*(S0.*PHI0-MEC.mb.PHIs(index));
          Rmbeff(index)=R0.*(1-S0);         
      end
      
      % process the non linear branches      
      for i=1:length(emnb)
          
         % get the branch number
         bn=emnb(i);
         pml=MEC.mb.pml{bn};
         nml=MEC.mb.nml{bn};
        
         % update system matix
         A(pml,pml)=A(pml,pml)+Rmbeff(bn);
         A(pml,nml)=A(pml,nml)-Rmbeff(bn);
         A(nml,pml)=A(nml,pml)-Rmbeff(bn);
         A(nml,nml)=A(nml,nml)+Rmbeff(bn);

         % update source vector
         b(pml)=b(pml)-Fsmbeff(bn);
         b(nml)=b(nml)+Fsmbeff(bn);
           
      end  % processing of nonlinear branches
 
      % solve the system
      if sum(sum(isinf(A)))||sum(sum(isnan(A)))
         k=MEC.mxit;
         x=zeros(size(b));
      else
         Nem=MEC.Nem;
         A11=A(1:Nem,1:Nem);
         A12=A(1:Nem,Nem+1:end);
         A21=A(Nem+1:end,1:Nem);
         A22=A(Nem+1:end,Nem+1:end);
         b1=b(1:Nem);
         b2=b(Nem+1:end);
         index=abs(A11)>0;
         NA11=mean(abs(A11(index)));
         index=abs(A22)>0;
         NA22=mean(abs(A22(index)));
         s=sqrt(NA11/NA22);
         sA=[A11 s*A12; s*A21 s*s*A22];
         sb=[b1; s*b2];
         x=sA\sb;
         x(Nem+1:end)=x(Nem+1:end)*s;
     end
      
     PHIm=x(1:MEC.Nem);
     Fn=x(MEC.Nem+1:end);
      
     % determine the branch fluxes
     for bn=1:MEC.Nemb
        pml=MEC.mb.pml{bn};
        nml=MEC.mb.nml{bn};
        PHImb(bn)=sum(PHIm(pml))-sum(PHIm(nml));
     end
  
     % determine the error
     e=norm(PHImht-PHIm)-(MEC.rec*0.5*norm(PHImht+PHIm)+MEC.aec);
     
     % update the iteration count
     k=k+1;

   end % while loop
 
   % test for convergence 
   c=(e<0)&&(k<=MEC.mxit);
  
end

% END Magnetically Non-Linear Problem-----------------------------------

% Solve for mesh branch information
if (MEC.Nemb>0)&&(nargout~=4)
   PHImR=PHImb-MEC.mb.PHIs;    
   Fmb=Fsmbeff+PHImb.*Rmbeff;
end

% Solve for nodal branch information
if (MEC.Nenb>0)&&(nargout~=5)
   Fnb=zeros(MEC.Nenb,1);
   pn=MEC.nb.pn;
   nn=MEC.nb.nn;
   index=pn>0;
   Fnb(index)=Fn(pn(index));
   index=nn>0;
   Fnb(index)=Fnb(index)-Fn(nn(index));
   PHInP=(Fnb-MEC.nb.Fs).*MEC.nb.P;
   PHInb=PHInP+MEC.nb.PHIs;
end

% output information
switch nargout
    case 4,
       varargout={Fn,Fnb,PHInb,c};
    case 5,
       varargout={PHIm,PHImb,PHImR,Fmb,c};
    case 9,
       varargout={PHIm,PHImb,PHImR,Fmb,Fn,Fnb,PHInb,PHInP,c};
    otherwise
       error('Inappropriate number of output arguments to mec_solve');
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