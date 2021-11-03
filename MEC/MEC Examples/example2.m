%---------------------------------------------------
%  SCRIPT:      example2.m
%  DESCRIPTION: A simple MEC of the EI core inductor
%               to demonstrate how to use the MEC
%               toolbox using mixed Mesh/Nodal Analysis
%  AUTHOR:      example 1 written by Cahya Harianto 
%               example 2 is example 1 modified by Scott D. Sudhoff                               
%               Purdue University
%               Electrical Engineering Building
%               465 Northwestern Avenue
%               West Lafayette, IN 47907-2035
%               sudhoff@ecn.purdue.edu
%               765-497-7648
%----------------------------------------------------

% Initialize
clear all; close all; clc

% Define some handy numbers
mm = 1e-3;
mu0 = 4*pi*10^-7;

% Specify the dimensions of the EI Core based on
% J.Cale, S.D.Sudhoff, and L.Tan, "Accurately Modeling EI Core Inductors 
% Using a High-Fidelity Magnetic Equivalent Circuit Approach", 
% IEEE Transactions on Magnetics, Vol.42, No.1, January 2006
we = 28.7 * mm;
wc = 57.4 * mm;
ws = 37.0 * mm;
ds = 48.8 * mm;
wb = 27.3 * mm;
wi = 30.5 * mm;
g  = 2.96 * mm;
d  = 119.8 * mm;
N  = 40;

% Curent range of interest
i  = linspace(0,100,100);
Npoints = length(i);

% Determine the MMF source on the circuit
F = N*i;

% Calculate the permeances on the air gap
P1 = mu0*(we*d)/g;
P2 = mu0*(wc*d)/g;
P3 = P1;

% 3C90 material parameters to calculate permeability as a function of B
M3C90.mur = 22340.9259;
M3C90.muB.a   = [  1.1542      0.049742    0.049644    0.041155];
M3C90.muB.b   = [431.1763       2.29503    15.04824    74.28908];
M3C90.muB.g   = [  0.4742        2.7955     0.59862     0.43996];
M3C90.muB.d   = M3C90.muB.a./M3C90.muB.b;
M3C90.muB.t   = exp(-M3C90.muB.b.*M3C90.muB.g);
M3C90.muB.h   = M3C90.muB.a.*M3C90.muB.t;
M3C90.muB.e   = M3C90.muB.t./(M3C90.muB.t+1);
M3C90.muB.z   = 1./(M3C90.muB.t+1);

% initial flux linkage and convergence
lambda=NaN(Npoints,1);
c=NaN(Npoints,1);

% Initialize the MEC
MEC = mec_init(4,10,5,3,1);

% Establish a list of non-linear materials
MEC = mec_nl_material(MEC,1,@muB,M3C90);

% Mesh Branches (MB)
MEC = mec_nls_branch(MEC,1,wc*d,(wb/2+ds/2),1,0,0.0,[1 -2],0.0); % MB 1                     
MEC = mec_nl_branch(MEC,2,we*d,(wb/2+ds/2),1,[1],0.0);           % MB 2      
MEC = mec_nl_branch(MEC,3,we*d,(wb/2+ds/2),1,[2],0.0);           % MB 3     
MEC = mec_nl_branch(MEC,4,wc*d,ds/2,1,[1 -2],0.0);               % MB 4
MEC = mec_nl_branch(MEC,5,we*d,ds/2,1,[1],0.0);                  % MB 5
MEC = mec_nl_branch(MEC,6,we*d,ds/2,1,[2],0.0);                  % MB 6
MEC = mec_nl_branch(MEC,7,wb*d,(wc/2+ws+we/2),1,[1],0.0);        % MB 7       
MEC = mec_nl_branch(MEC,8,wb*d,(wc/2+ws+we/2),1,[2],0.0);        % MB 8
MEC = mec_nl_branch(MEC,9,wi*d,(wc/2+ws+we/2),1,[3],0.0);        % MB 9
MEC = mec_nl_branch(MEC,10,wi*d,(wc/2+ws+we/2),1,[4],0.0);       % MB 10
                          
% Explicit Nodes (EN)
MEC = mec_enode(MEC,1,[1]);                                      % EN 1
MEC = mec_enode(MEC,2,[-3]);                                     % EN 2
MEC = mec_enode(MEC,3,[3 -4]);                                   % EN 3
MEC = mec_enode(MEC,4,[4]);                                      % EN 4
MEC = mec_enode(MEC,5,[-2]);                                     % EN 5
    
% Nodal Branches (NB)
MEC = mec_lnb(MEC,1,1,2,P1,0,0);                                 % NB 1
MEC = mec_lnb(MEC,2,3,0,P2,0,0);                                 % NB 2
MEC = mec_lnb(MEC,3,4,5,P3,0,0);                                 % NB 3

% Test contruction of MEC
if (mec_test(MEC)~=1)
   error('Bad MEC');
end

% Loop over current values
for n = 1:Npoints  

   % Assign MMF
   MEC.mb.Fs(1)=F(n);   
   
   % Solve the MEC
   [PHIm,PHImb,PHImR,Fmb,Fn,Fnb,PHInb,PHInP,c(n)] = mec_solve(MEC);

   % Compute the flux linkage
   lambda(n) = - PHImR(1)*N;
       
end

% Generate the lambda-vs-current curve
if (min(c)<1)
   error('Solution Did Not Converge');
end

% Plot flux linkage versus current
plot(i,lambda)
xlabel('i, (A)')
ylabel('\lambda, (Vs)')
grid on
axis([0 100 0 0.2]);

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
%  License along with MEC Toolbox. If not, 
%  see <http://www.gnu.org/licenses/>.