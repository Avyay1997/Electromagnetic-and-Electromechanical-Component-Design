%  FILE:        awg_floor.m
%  FUNCTION:    [WP]=awg_floor(rw)
%  DESCRIPTION: Assigns wire parameters so that it is the closest AWG with
%               a radius less than the corresponding AWG parameters              
%
%  INPUTS:      rw        - desired radius
%
%  OUTPUTS:     WP.desc   - conductor desription
%               WP.area   - wire area (m^2)
%               WP.diam   - wire diameter (m)
%               WP.rad    - wire radius (m)
%
%  REFERENCE:   Herrington, D. E., Handbook of Electronic Tables and
%               Formulas. Indianapolis, IN: Howard W. Sams & Co., 1959.

function [WP]=awg_floor(rw)

%  Define AWG diameter
diam=[...
    11.68400 ...
    10.40384 ...
     9.26592 ...
     8.25246 ...
     7.34822 ...
     6.54304 ...
     5.82676 ...
     5.18922 ...
     4.62026 ...
     4.11480 ...
     3.66522 ...
     3.26390 ...
     2.90576 ...
     2.58826 ...
     2.30378 ...
     2.05232 ...
     1.82880 ...
     1.62814 ...
     1.45034 ...
     1.29032 ...
     1.15062 ...
     1.02362 ...
     0.91186 ...
     0.81280 ...
     0.72390 ...
     0.64516 ...
     0.57404 ...
     0.51054 ...
     0.45466 ...
     0.40386 ...
     0.36068 ...
     0.32004 ...
     0.28702 ...
     0.25400 ...
     0.22606 ...
     0.20320 ...
     0.18034 ...
     0.16002 ...
     0.14224 ...
     0.12700 ...
     0.11430 ...
     0.10160 ...
     0.08890 ...
     0.07874]*1e-3;

%  Define AWG gauge
gauge=['4/0';'3/0';'2/0';'1/0';' 1 ';' 2 ';' 3 ';' 4 ';' 5 '; ...
       ' 6 ';' 7 ';' 8 ';' 9 ';' 10';' 11';' 12';' 13';' 14'; ...
       ' 15';' 16';' 17';' 18';' 19';' 20';' 21';' 22';' 23'; ...
       ' 24';' 25';' 26';' 27';' 28';' 29';' 30';' 31';' 32'; ...
       ' 33';' 34';' 35';' 36';' 37';' 38';' 39';' 40'];

%  Find AWGs that are less than the desired wire radius   
index=find(2*rw<diam);
diam=diam(index);
gauge=gauge(index,:);

%  Define custom wire parameters if desired wire radius is greater than AWG
if isempty(index)
   WP.desc='CUS';
   WP.diam=2*rw;
   WP.rad=rw;
   WP.area=pi*rw^2;
   
%  Determine AWG closest to desired wire diameter   
else 
   error=abs(diam.^2 - (2*rw)^2);
   [minerror,index]=min(error);
   
   %  Define outputs
   WP.desc=gauge(index,:);
   WP.diam=diam(index);
   WP.rad=WP.diam/2.0;
   WP.area=pi*(WP.rad).^2;
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