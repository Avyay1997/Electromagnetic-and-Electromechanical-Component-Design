function [SP]=steel_catalog(n)
% steel_catalag assigns steel parameters to a parameter structure
% 
% Call:
% SP=steel_catalog(n)
%
% Inputs:
% n              = steel number
%                  n=1 - M19 Steel
%                  n=2 - M36 Steel
%                  n=3 - M43 Steel
%                  n=4 - M47 Steel
%                  n=5 - Hiperco 50 Steel
%
% Outputs:
% SP             = structure of ferrite parameters 
%  SP.desc       = material description
%  SP.Msat       = saturated magnetization (T)
%  SP.mur        = relative permeability in linear region
%  SP.Blim       = recommended limit on B to avoid saturation (T) 
%                 (this is taken to be the point where absolute relative 
%                  permeability hits 1000)
%  SP.row        = mass density (kg/m^3)
%  SP.BH         = structure of BH curve parameters
%   SP.BH.m      = magnetization coefficients (T)
%   SP.BH.n      = exponents
%   SP.BH.h      = field intensity breakpoints (A/m)
%  SP.MSE        = structure of MSE loss parameters
%   SP.MSE.alpha = frequency exponent
%   SP.MSE.beta  = flux density exponent
%   SP.MSE.kh    = hysteresis loss coefficient (J/m^3)
%   SP.MSE.ke    = eddy current loss coefficient (J*s/m^3)
%  SP.muB        = structure of mu(B) parameters
%   SP.muB.mur   = initial relative permeability of anhysteretic curve
%   SP.muB.a     = vector of alpha coefficients (1/T)
%   SP.muB.b     = vector of beta exponential coefficients (1/T)
%   SP.muB.g     = vector of gamma exponential offsets (T)
%   SP.muB.d     = vector of delta coefficients
%   SP.muB.e     = vector of epsilon values
%   SP.muB.z     = vector of zeta values
%   SP.muB.h     = vector of eta values (1/T)
%   SP.muB.t     = vector of theta values 
%
% Written by:
% Grant Shane for Scott D. Sudhoff                               
% School of Electrical and Computer Engineering
% 1285 Electrical Engineering Building
% West Lafayette, IN 47907-1285
% Phone: 765-494-3246
% Fax: 765-494-0676
% E-mail: sudhoff@ecn.purdue.edu
%
% References:
% ASM International Metals Handbook, Properties and Selection:
% Nonferrous Alloys and Special-Purpose Materials.  vol. 2, 
% 10th ed., 1990.
%
% Mapes & Sprowl Steel, Elk Grove Village, IL, private communication, 2010.
%
% Shane, G. and S. D. Sudhoff,  "Refinements in Anhysteretic
% Characterization and Permeability Modeling,"  IEEE
% Transactions on Magnetics, vol. 46, no. 11 November 2010.
%
% Modifications:
% This routine was modified on 12/30/2010 in order to fix an error.  In
% particular, the turns ratio was not incorporated into the loss
% calculation of the MSE Parameters.  This routine contains the updated 
% MSE Parameters.
%
% Note: This code was taken form the Power Magnetic Material Toolbox by
%       S.D. Sudhoff

switch n
    
    case 1 % M19-----------------------------------------------------------
        
        % general parameters
        SP.desc='M19';
        SP.Msat=1.4311;
        SP.mur=32565.7512;
        SP.Blim=1.3922;
        SP.row=7402;
        SP.k=16.7; % [ASM, pg. 767]
        SP.c=469; % [Mapes & Sprowl]
        
        % BH parameters
        SP.BH.m = [1.4311    -0.18394     0.15175    0.051208];
        SP.BH.n = [1      3.9385           5      3.3648];
        SP.BH.h = [34.14363      85.25885      151.4965      311.3874];
        
        % MSE parameters
        SP.MSE.alpha=1.3375;
        SP.MSE.beta=1.8167;
        SP.MSE.kh=50.691;
        SP.MSE.ke=2.7513e-2;
        
        % mu(B) parameters
        SP.muB.mur = 32685.6784;
        SP.muB.a   = [0.098611   0.0014823    0.001435    0.001435];
        SP.muB.b   = [69.73973      1.949541      162.2767      3.598553];
        SP.muB.g   = [1.399      2.1619      1.2475      2.0377];
        
    case 2 % M36-----------------------------------------------------------
        
        % general parameters
        SP.desc='M36';
        SP.Msat=1.3652;
        SP.mur=27402.6656;
        SP.Blim=1.3371;
        SP.row=7018;
        SP.k=18.8; % [ASM, pg. 767]
        SP.c=465; % [Mapes & Sprowl]
        
        % BH parameters
        SP.BH.m = [1.3652    -0.40925     0.22662    -0.28663];
        SP.BH.n = [1      2.0941      2.5541      4.8654];
        SP.BH.h = [24.47644      20.38454      110.6466      86.54568];
        
        % MSE parameters
        SP.MSE.alpha=1.3403;
        SP.MSE.beta=1.7991;
        SP.MSE.kh=64.064;
        SP.MSE.ke=4.0473e-2;
        
        % mu(B) parameters
        SP.muB.mur = 26673.3745;
        SP.muB.a   = [0.22599    0.043195    0.031118   0.0043748];
        SP.muB.b   = [271.8443      97.31738      42.29465     0.8058081];
        SP.muB.g   = [1.35065           10      1.32415      5.38168]; 
        
    case 3 % M43-----------------------------------------------------------
        
        % general parameters
        SP.desc='M43';
        SP.Msat=1.5977;
        SP.mur=24817.2695;
        SP.Blim=1.3902;
        SP.row=7291;
        SP.k=20.9; % [ASM, pg. 767]
        SP.c=461; % [Mapes & Sprowl]
        
        % BH parameters
        SP.BH.m = [1.5977    -0.19905    -0.27081     0.30728];
        SP.BH.n = [1      4.8746       1.419      2.9265];
        SP.BH.h = [50.00814      64.59417      843.8071      116.3637];
        
        % MSE parameters
        SP.MSE.alpha=1.2785;
        SP.MSE.beta=1.7543;
        SP.MSE.kh=85.036;
        SP.MSE.ke=4.1409e-2;
        
        % mu(B) parameters
        SP.muB.mur = 24891.9232;
        SP.muB.a   = [0.072885   0.0039565   0.0025477       0.001];
        SP.muB.b   = [33.9743      1.14108      3.37941      48.9292];
        SP.muB.g   = [1.4138      4.2848      9.9998      1.4731];
        
    case 4 % M47-----------------------------------------------------------
        
        % general parameters
        SP.desc='M47';
        SP.Msat=1.8571;
        SP.mur=9846.0807;
        SP.Blim=1.4874;
        SP.row=7585;
        SP.k=37.7; % [ASM, pg. 767]
        SP.c=456; % [Mapes & Sprowl]
        
        % BH parameters
        SP.BH.m = [1.8571     0.83433    -0.51725      1.1778];
        SP.BH.n = [1       3.304      1.2606      2.0209];
        SP.BH.h = [289.999       159.604      160.2678      296.7384];
        
        % MSE parameters
        SP.MSE.alpha=1.2559;
        SP.MSE.beta =1.685;
        SP.MSE.kh=149.0;
        SP.MSE.ke=0.261;
        
        % mu(B) parameters
        SP.muB.mur = 9875.0004;
        SP.muB.a   = [0.050919     0.03933       0.001       0.001];
        SP.muB.b   = [18.03839      20.61242     0.9336389      115.8623];
        SP.muB.g   = [1.5613      3.3712      5.2715      1.3996];
        
    case 5 % Hiperco 50----------------------------------------------------
        
        % general parameters
        SP.desc='Hiperco50';
        SP.Msat=2.4431;
        SP.mur=42891.5684;
        SP.Blim=2.0709;
        SP.row=7844.64; % [measured, Datasheet: 8110 kg/m/K]
        SP.k=29.8064;   % [Carpenter Technical Datasheet]
        SP.c=418.4;     % [phone conversation: 0.1 Btu/lb/F]
        
        % BH parameters
        SP.BH.m = [2.4431    -0.43856     0.19526    -0.82542];
        SP.BH.n = [1      1.1011      4.9986      2.2891];
        SP.BH.h = [34.85781      54.62708      253.8951      92.44246];
        
        % MSE parameters
        SP.MSE.alpha=1.0839;
        SP.MSE.beta =1.855;
        SP.MSE.kh=74.528;
        SP.MSE.ke=3.0185e-2;
        
        % mu(B) parameters
        SP.muB.mur = 43371.9609;
        SP.muB.a   = [0.43708   0.0003068  0.00026279  0.00024516];
        SP.muB.b   = [17.13367      2.139356      163.4348      1.476588];
        SP.muB.g   = [2.2836      1.3692      1.6772       3.494];
        
    otherwise %------------------------------------------------------------
        
        % assign NaN to all parameters and exit function
        SP.desc='BAD';
        SP.Msat=NaN;
        SP.mur=NaN;
        SP.Blim=NaN;
        SP.row=NaN;
        SP.k=NaN;
        SP.c=NaN;
        
        SP.BH.m=NaN;
        SP.BH.h=NaN;
        SP.BH.n=NaN;
        
        SP.MSE.alpha=NaN;
        SP.MSE.beta=NaN;
        SP.MSE.kh=NaN;
        SP.MSE.ke=NaN;
        
        SP.muB.mur=NaN;
        SP.muB.a=NaN;
        SP.muB.b=NaN;
        SP.muB.g=NaN;
        SP.muB.d=NaN;
        SP.muB.t=NaN;
        SP.muB.h=NaN;
        SP.muB.e=NaN;
        SP.muB.z=NaN;
        return     
        
end

% calculated mu(B) parameters
SP.muB.d=SP.muB.a./SP.muB.b;
SP.muB.t=exp(-SP.muB.b.*SP.muB.g);
SP.muB.h=SP.muB.a.*SP.muB.t;
SP.muB.e=SP.muB.t./(SP.muB.t+1);
SP.muB.z=1./(SP.muB.t+1);

end

%  Copyright 2013 - Scott Sudhoff 
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