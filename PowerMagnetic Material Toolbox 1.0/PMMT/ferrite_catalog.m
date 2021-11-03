%  FILE:        ferrite_catalog.m
%  FUNCTION:    [FP]=ferrite_catalog(n)
%  DESCRIPTION: This routine assigns ferrite parameters.
%
%  INPUTS:      n         - ferrite number
%               n=1       - MN8CX
%               n=2       - MN60LL
%               n=3       - MN67
%               n=4       - MN80C
%               n=5       - 3C90
%
%  OUTPUTS:     FP        - structure of ferrite parameters 
%                  FP.desc      - material description
%                  FP.Msat      - saturated magnetization (T)
%                  FP.mur       - relative permeability in linear region
%                                 (from anhysteretic results)
%                  FP.Blim      - recommended limit on B to avoid saturation (T) 
%                                 (This is point where absolute relative 
%                                 permeability hits 1000)
%                  FP.row       - mass density (kg/m^3)
%
%               FP.BH     - structure of BH curve parameters
%                  FP.BH.m      - magnetization coefficients (T)
%                  FP.BH.n      - exponents
%                  FP.BH.h      - field intensity breakpoints (A/m)
%
%               FP.MSE    - structure of MSE loss parameters
%                  FP.MSE.alpha - frequency exponent
%                  FP.MSE.beta  - flux density exponent
%                  FP.MSE.kh    - hysteresis loss coefficient (J/m^3)
%                  FP.MSE.ke    - eddy current loss coefficient (J*s/m^3)
%
%               FP.muB    - structure of mu(B) parameters
%                  FP.muB.mur   - initial relative permeability of
%                                 anhysteretic curve
%                  FP.muB.a     - vector of alpha coefficients (1/T)
%                  FP.muB.b     - vector of beta exponential coefficients (1/T)
%                  FP.muB.g     - vector of gamma exponential offsets (T)
%                  FP.muB.d     - vector of delta coefficients
%                  FP.muB.e     - vector of epsilon values
%                  FP.muB.z     - vector of zeta values
%                  FP.muB.h     - vector of eta values (1/T)
%                  FP.muB.t     - vector of theta values 
%
%  AUTHOR:      Grant Shane for Scott D. Sudhoff                               
%               Purdue University
%               Electrical Engineering Building
%               465 Northwestern Avenue
%               West Lafayette, IN 47907-2035
%               sudhoff@ecn.purdue.edu
%               765-497-7648
%
%  REFERENCE:  Shane, G. and S. D. Sudhoff,  "Refinements in Anhysteretic
%              Characterization and Permeability Modeling,"  IEEE
%              Transactions on Magnetics, submitted for publication.

function [FP]=ferrite_catalog(n)

switch n
    
    case 1 % MN8CX---------------------------------------------------------
        
        % general parameters
        FP.desc='MN8CX';
        FP.Msat=0.5996;
        FP.mur=6376.9879;
        FP.Blim=0.41782;
        FP.row=4611.7;
        
        % BH parameters
        FP.BH.m = [0.5996    -0.13838    -0.07605     0.26724];
        FP.BH.n = [1      1.1001      1.3202      2.1278];
        FP.BH.h = [109.41615       293.1435      1312.9303      87.234079];
        
        % MSE parameters
        FP.MSE.alpha=1.01;
        FP.MSE.beta=2.70;
        FP.MSE.kh=21.8;
        FP.MSE.ke=0;
        
        % mu(B) parameters
        FP.muB.mur = 6325.9986;
        FP.muB.a   = [0.41961     0.18748    0.010849   0.0081058];
        FP.muB.b   = [269.1083       983.769       8.69827      155.8723];
        FP.muB.g   = [0.42322     0.45413     0.60172     0.38529];
        
    case 2 % MN60LL--------------------------------------------------------
   
        % general parameters
        FP.desc='MN60LL';
        FP.Msat=0.456;
        FP.mur=15761.7796;
        FP.Blim=0.44428;
        FP.row=4818.5;
        
        % BH parameters
        FP.BH.m = [0.456   -0.072531     0.1192   -0.12438];
        FP.BH.n = [1      4.6134      1.7099      1.4343];
        FP.BH.h = [10.9362  24.345  33.3318  5.53123];
        
        % MSE parameters
        FP.MSE.alpha=1.0342;
        FP.MSE.beta=2.301;
        FP.MSE.kh=6.7212;
        FP.MSE.ke=1.4265e-12;

        % mu(B) parameters
        FP.muB.mur = 15042.4481;
        FP.muB.a   = [2.1948    0.021436   0.0078737    0.001005];
        FP.muB.b   = [306.0455      998.3986      180.4701      6.824599];
        FP.muB.g   = [0.453      8.2682     0.40563     0.38732];
        
    case 3 % MN67----------------------------------------------------------
        
        FP.desc='MN67';
        FP.Msat=0.57499;
        FP.mur=26999.419;
        FP.Blim=0.46327;
        FP.row=4794.7;
        
        % BH parameters
        FP.BH.m = [0.57499      0.4721     0.24465    -0.11526];
        FP.BH.n = [1      1.7699      1.7503      1.3535];
        FP.BH.h = [149.4072      16.83932      81.27764      119.3392];
        
        % MSE parameters
        FP.MSE.alpha=1.035;
        FP.MSE.beta=1.667;
        FP.MSE.kh=11.2;
        FP.MSE.ke=9.93e-3;
        
        % mu(B) parameters
        FP.muB.mur = 27084.7572;
        FP.muB.a   = [1.0588    0.021505    0.021505       0.001];
        FP.muB.b   = [74.31812      9.633883      48.52124      106.7594];
        FP.muB.g   = [0.55342      0.7891      0.4303      0.3785];
        
    case 4 % MN80C---------------------------------------------------------
        
        % general parameters
        FP.desc='MN80C';
        FP.Msat=0.70799;
        FP.mur=12399.9358;
        FP.Blim=0.44343;
        FP.row=4709.8;
        
        % BH parameters
        FP.BH.m = [0.70799    -0.47286    -0.12588     0.38982];
        FP.BH.n = [1      1.1406      2.1121      1.4594];
        FP.BH.h = [39.77445      76.97549      183.2491      84.52192];
        
        % MSE parameters
        FP.MSE.alpha=1.002;
        FP.MSE.beta=2.114;
        FP.MSE.kh=27.72;
        FP.MSE.ke=0;
        
        % mu(B) parameters
        FP.muB.mur = 12402.0239;
        FP.muB.a   = [0.76926    0.048976     0.01336   0.0028669];
        FP.muB.b   = [85.6801      16.7675      97.5145      2.98234];
        FP.muB.g   = [0.49503     0.58504     0.41129      1.0098];
        
    case 5 % 3C90----------------------------------------------------------
        
        % general parameters
        FP.desc='3C90';
        FP.Msat=0.50163;
        FP.mur=22340.9259;
        FP.Blim=0.44916;
        FP.row=4743.3;
        
        % BH parameters
        FP.BH.m = [0.50163     0.02493   -0.060579    0.018615];
        FP.BH.n = [1      1.9128      1.2282      2.8616];
        FP.BH.h = [17.987023      29.251923      88.639093      1080.3857];
        
        % MSE parameters
        FP.MSE.alpha=1.0064;
        FP.MSE.beta=1.9272;
        FP.MSE.kh=84.715;
        FP.MSE.ke=.072189e-3;
        
        % mu(B) parameters
        FP.muB.mur = 22266.6178;
        FP.muB.a   = [1.1542    0.049742    0.049644    0.041155];
        FP.muB.b   = [431.1763       2.29503      15.04824      74.28908];
        FP.muB.g   = [0.4742      2.7955     0.59862     0.43996];
        
    otherwise %------------------------------------------------------------
        
        % assign NaN to all parameters and exit function
        FP.desc='BAD';
        FP.Msat=NaN;
        FP.mur=NaN;
        FP.Blim=NaN;
        FP.row=NaN;
        
        FP.BH.m=NaN;
        FP.BH.h=NaN;
        FP.BH.n=NaN;
        
        FP.MSE.alpha=NaN;
        FP.MSE.beta=NaN;
        FP.MSE.kh=NaN;
        FP.MSE.ke=NaN;
        
        FP.muB.mur=NaN;
        FP.muB.a=NaN;
        FP.muB.b=NaN;
        FP.muB.g=NaN;
        FP.muB.d=NaN;
        FP.muB.t=NaN;
        FP.muB.h=NaN;
        FP.muB.e=NaN;
        FP.muB.z=NaN;
        return 
        
end

% calculated mu(B) parameters
FP.muB.d=FP.muB.a./FP.muB.b;
FP.muB.t=exp(-FP.muB.b.*FP.muB.g);
FP.muB.h=FP.muB.a.*FP.muB.t;
FP.muB.e=FP.muB.t./(FP.muB.t+1);
FP.muB.z=1./(FP.muB.t+1);

end

%  Copyright 2010 - Scott Sudhoff 
% 
%  The program is distributed under the terms of the GNU Lesser
%  General Public License(GNU LGPL). 
% 
%  This file is part of the Power Magnetic Material Toolbox.
% 
%  The Power Magnetic Material Toolbox is free software: 
%  you can redistribute it and/or modify
%  it under the terms of the GNU Lesser General Public License as published
%  by the Free Software Foundation, either version 3 of the License, or 
%  (at your option) any later version.
% 
%  The Power Magnetic Material Toolbox is distributed 
%  in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
%  See the GNU Lesser General Public License for more details.
% 
%  You should have received a copy of the GNU Lesser General Public
%  License along with the Power Magnetic Material Toolbox.  If not, see
%  <http://www.gnu.org/licenses/>.