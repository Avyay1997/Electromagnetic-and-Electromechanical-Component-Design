%---------------------------------------------------
%  SCRIPT:      example1.m
%  DESCRIPTION: A simple MEC of the EI core inductor
%               to demonstrate how to use the MEC
%               toolbox
%  AUTHOR:      Cahya Harianto for Scott D. Sudhoff                               
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
F1 = N*i;

% Calculate the permeances on the air gap
R11 = g/(mu0*wc*d);
R12 = g/(mu0*we*d);
R13 = R12;

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

% Initial flux linkage and convergence
lambda=NaN(Npoints,1);
c=NaN(Npoints,1);

% Initialize the MEC
MEC = mec_init(2,13,1);    
   
% Establish a list of non-linear materials
MEC = mec_nl_material(MEC,1,@muB,M3C90);

% Branches
MEC= mec_nls_branch(MEC,1,wc*d,(wb/2+ds/2),1,0.0,0.0,[1 -2],0.0);  % MB 1
MEC= mec_nl_branch(MEC,2,we*d,(wb/2+ds/2),1,[1],0.0);              % MB 2
MEC= mec_nl_branch(MEC,3,we*d,(wb/2+ds/2),1,[2],0.0);              % MB 3
MEC= mec_nl_branch(MEC,4,wc*d,ds/2,1,[1 -2],0.0);                  % MB 4
MEC= mec_nl_branch(MEC,5,we*d,ds/2,1,[1],0.0);                     % MB 5
MEC= mec_nl_branch(MEC,6,we*d,ds/2,1,[2],0.0);                     % MB 6
MEC= mec_nl_branch(MEC,7,wb*d,(wc/2+ws+we/2),1,[1],0.0);           % MB 7
MEC= mec_nl_branch(MEC,8,wb*d,(wc/2+ws+we/2),1,[2],0.0);           % MB 8
MEC= mec_nl_branch(MEC,9,wi*d,(wc/2+ws+we/2),1,[1],0.0);           % MB 9 
MEC= mec_nl_branch(MEC,10,wi*d,(wc/2+ws+we/2),1,[2],0.0);          % MB10
MEC= mec_lr_branch(MEC,11,R11,[1 -2]);                             % MB11
MEC= mec_lr_branch(MEC,12,R12,[1]);                                % MB12
MEC= mec_lr_branch(MEC,13,R13,[2]);                                % MB13

% Loop over current values
for n = 1:Npoints

   % Assign MMF
   MEC.mb.Fs(1)=F1(n);                   

   % Solve the MEC
   [Pm,Pb,PR,Fb,converge] = mec_solve(MEC);
   
   % Compute total flux linkage
   lambda(n) = -PR(1)*N;
   
   % Copy convergence
   c(n)=converge;
    
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
%  License along with MEC Toolbox.  
%  If not, see <http://www.gnu.org/licenses/>.
