% pmac_design is a multi-objective optimal design of a permanent magnet
%             ac machine. Used in Chapter 9 Section 11 of "Power Magnetic 
%             Devices: A Multi-Objective Design Approach" by S.D. Sudhoff
%
% Version/Date:
% June 12, 2013
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
% Email: asah@purdue.edu

% read in common routines--------------------------------------------------
path(path,'Common Routines');

% read design specifications-----------------------------------------------
pmac_specs;

% determine parameters for genetic algorithm-------------------------------
Nobj=2;                              % Nobj = number of objectives
Obj=0;                               % Obj = objective to optimize
                                     %   1 = minimize size
                                     %   2 = minimize loss
                                     %   0 = multiobjective optimization
Npop=3000;                           % Population Size
Ngen=3000;                           % Number of generations

GAP=gapdefault(Nobj,Obj,Npop,Ngen);
GAP.op_list=[1];
GAP.ev_pp=true;
GAP.rp_lvl=1;
GAP.mg_nreg=round(GAP.fp_npop/400);

% gene description---------------------------------------------------------
%
%          min    max  chrm  chrm      par
%          val    val  type   id         #     description
GAP.gd =[   1       4    1    1;  ... %  1  stator steel type
            1       4    1    1;  ... %  2  rotor steel type
            1       2    1    1;  ... %  3  conductor type [MOD]
            1       7    1    1;  ... %  4  permanent magnet type
            4       8    1    1;  ... %  5  pole pairs
            0     10*cm  2    1;  ... %  6  depth of inert region(m)
           1*mm    5*cm  3    1;  ... %  7  depth of rotor backiron (m)
           1*mm    5*cm  3    1;  ... %  8  depth of pm (m)
          0.5*mm   2*mm  2    1;  ... %  9  air gap 
           1*mm    5*cm  3    1;  ... % 10  depth of tooth base (m)
          0.05     0.95  2    1;  ... % 11  tooth fraction
           1*mm    5*cm  3    1;  ... % 12  depth of stator backiron(m)
          0.05     0.95  2    1;  ... % 13  pm fraction 
           1*cm   50*cm  3    1;  ... % 14  active length
           10      1e3   3    1;  ... % 15  pk phs cond.density(cond/rad)
           0.1     0.7   2    1;  ... % 16  relative 3rd harmon cond.
                                  ... %     density
            1     1.15    3    1; ... % 17  Kt1: torque multiplier
            1     1.15    3    1; ... % 18  Kt2: torque multiplier
            1     1.15    3    1; ... % 19  Kt3: torque multiplier
            1     1.15    3    1; ... % 20  Kt4: torque multiplier
            1     1.15    3    1; ... % 21  Kt5: torque multiplier
            1     1.15    3    1]; ... % 22 Kt6: torque multiplier
           
     
% conduct the optimization-------------------------------------------------
[fP,GAS,bestparameters]=gaoptimize(@pmac_fitness,GAP,D);
save resultsPMAC3000