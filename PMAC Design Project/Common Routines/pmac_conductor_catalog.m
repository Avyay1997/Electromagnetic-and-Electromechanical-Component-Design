function [C]=pmac_conductor_catalog(n)
% pmac_conductor_catalog creates of structure of conductor data.
%                        "Power Magnetic Devices: A Multi-Objective Design 
%                        Approach," by S.D. Sudhoff
%
% Call:
% C=pmac_conductor_catalog(n)
%
% Inputs:
% n        = conductor code
%           n=1       - copper
%           n=2       - aluminum
%
% Outputs:
% C        = structure of conductor parameters
%  C.desc  = conductor description
%  C.sigma0= conductivity at t0 (1/(Ohm-m))
%  C.sigmac= same as CP.sigma0 (backwards compatibility)
%  C.Jlim  = recommended maximum current density (A/m^2) 
%  C.row   = mass density (kg/m^3)
%
% Written for:
% S.D. Sudhoff
% School of Electrical and Computer Engineering
% 1285 Electrical Engineering Building
% West Lafayette, IN 47907-1285
% Phone: 765-494-3246
% Fax: 765-494-0676
% E-mail: sudhoff@ecn.purdue.edu
%
% References:
%
% Fink, D. G. (ed.), and H. W. Beaty (ed.), Standard Handbook 
% for Electrical Engineers. 13th ed. New York: McGraw-Hill, Inc., 1993.
%
% Lide, D. E. (ed.), CRC Handbook of Chemistry and Physics.
% 84th ed. Boca Raton, FL: CRC Press, 2003.
%
% Serway, R. A., Principles of Physics. 2nd ed. Fort Worth,
% TX:  Saunders College Pub., 1998.
%
% Stauffer, H. B., Engineer's Guide to the National Electric
% Code. Boston, MA: Jones and Bartlett Publishers, Inc., 2008.
%
% Young, H. D., University Physics.  7th ed. Reading, MA: 
% Addison-Wesley, 1992.
%
% Modified by Avyay Sah
% Email: asah@purdue.edu


switch n
   case 1
      C.desc='Copper';
      C.sigma0=3.94466874e7; % sigma for 150C  [Lide, sec. 12]
      C.Jlim=7.6025e6;  % calculated using 10 AWG at 90 C (194 F) which 
                        % has a max ampacity of 40 A [Stauffer, Table 3.3]
      C.row=8890;       % @ 20 C (pure copper, rolled, forged, or drawn 
                        % and then annealed) [Fink, pg. 4-3]
   case 2
      C.desc='Aluminum';
      C.sigma0=2.493056474e7; % sigma for 150C   [Lide, sec. 12]
      C.Jlim=6.6522e6;  % calculated using 10 AWG at 90 C (194 F) which 
                        % has a max ampacity of 35 A [Stauffer, Table 3.3]
      C.row=2705;       % @ 20 C (commercially hard-drawn) [Fink, pg. 4-3]
   otherwise
      error('Unknown conductor type');    
              
end % case
         
end