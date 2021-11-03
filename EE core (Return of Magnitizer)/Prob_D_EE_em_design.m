%  Prob_D_EE_em_design is a multi objective optimal design of a EE core
%               electromagnet based on an MEC analysis. 
%
% Version/Date:
% May 22, 2012
%
% Written by:
% S.D. Sudhoff                               
% Purdue University
% Electrical Engineering Building
% 465 Northwestern Avenue
% West Lafayette, IN 47907-2035
% sudhoff@ecn.purdue.edu
% 765-497-7648
%
% Modified by Avyay Sah
% Email asah@purdue.edu

% genetic algorithm parameters
ngen=2000;
npop=2000;
GAP=gapdefault(2,0,npop,ngen);
GAP.dv_act=0;
%GAP.ev_pp=true;
GAP.rp_lvl=0;
GAP.mg_nreg=10;

% design specifications
D.g=5e-3;                       % air gap (m)
D.Bmin=1.75;                    % flux density (T)
D.SZmax=10;                     % maximum allowed size sqrt(kg m^3)
D.PLmax=500;                    % maximum allowed loss (W)
D.kpf=0.6;                      % packing factor
D.cgd=5e-3;                     % vertical air gap clearance(m)
D.lc=50e-3;                     % length of core (m)
D.wc=50e-3;                     % width of core (m)

% gene description---------------------------------------------------------
%
%         min   max chrm  chrm      par
%         val   val type   id         #     description
GAP.gd= [2e-3  9e-1   3    1;  ...  %  1   wb  % width of base of core
         0.5    1.5   2    1;  ...  %  2   re  % ratio of width of leg of core
         2e-3  3e-1   3    1;  ...  %  3   ww  % width of conductor
         2e-3  3e-1   3    1;  ...  %  4   dw  % depth of conductor
          1    5e6    3    1;  ...  %  5   J   % current density
         5e-3  1e-1   3    1;  ...  %  6   cr  % corner clearnace
         5e-3  1e-1   3    1;  ...  %  7   cso % outside core clearance
         5e-3  1e-1   3    1];  ... %  8   cb  % base core clearance
          
[fP,GAS,bi,bf]= gaoptimize(@Prob_D_EE_em_fit,GAP,D);
 
save Prob_D_EE_em_design_results_2000