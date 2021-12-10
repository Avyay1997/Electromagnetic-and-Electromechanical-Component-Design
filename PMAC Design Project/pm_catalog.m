%  FILE:        pm_catalog.m
%  FUNCTION:    [MP]=pm_catalog(n)
%  DESCRIPTION: This routine assigns permanent magnet parameters.
%
%  INPUTS:      n         - permanent magnet number
%               n=1       - NdFeB-4SB
%               n=2       - NdFeB-30
%               n=3       - NdFeB-48
%               n=4       - SmCo-B15S
%
%  OUTPUTS:     MP        - data structure of permanent magnet parameters
%                  MP.desc   - magnet desription
%                  MP.brm    - residual flux density (T)
%                  MP.xpm    - susceptibility 
%                  MP.hmmn   - minimum field intensity before demag (A/m)
%                  MP.row    - mass density (kg/m^3)
%
%  AUTHOR:      Scott D. Sudhoff                               
%               Purdue University
%               Electrical Engineering Building
%               465 Northwestern Avenue
%               West Lafayette, IN 47907-2035
%               sudhoff@ecn.purdue.edu
%               765-497-7648

function [MP]=pm_catalog(n)

switch n
   case 1
      MP.desc='NdFeB-4SB';
      MP.brm=0.367;
      MP.xpm=0.2;
      MP.hmmn=-122e3;
      MP.row=6005;
   case 2
      MP.desc='NdFeB-30';
      MP.brm=1.1;
      MP.xpm=0.0;
      MP.hmmn=-436e3;
      MP.row=7389;
   case 3
      MP.desc='NdFeB-48';
      MP.brm=1.45;
      MP.xpm=0.1;
      MP.hmmn=-526e3;
      MP.row=7500;
   case 4
      MP.desc='SmCo-B15S';
      MP.brm=0.785;
      MP.xpm=0.1;
      MP.hmmn=-284e3;
      MP.row=7000;
   otherwise
      MP.desc='BAD';
      MP.brm=NaN;
      MP.xpm=NaN;
      MP.hmmn=NaN;
      MP.row=NaN;
end
          
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