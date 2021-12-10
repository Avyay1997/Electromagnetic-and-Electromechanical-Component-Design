% pmac_specs defines the design specifications for a multi-objective 
%            optimal design of a permanent magnet ac machine. Used in 
%            Chapter 9 Section 11 of "Power Magnetic Devices: A 
%            Multi-Objective Design Approach" by S.D. Sudhoff
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

% initialize---------------------------------------------------------------
close all
clear all

% handy numbers------------------------------------------------------------
rpm_to_radps=2*pi/60;                % convert RPM to rad/s
cm=1e-2;                             % centimeter in meters
mm=1e-3;                             % millimeter in meters

% Given Speeds for 6 operating points
a_cpsr = 3;                            % speed ratio
wrm_max = 20000*rpm_to_radps;
wrm_crn = wrm_max/a_cpsr;
wrm_nom = sqrt(wrm_crn*wrm_max);
WRM = [wrm_crn wrm_crn wrm_max wrm_crn/2 wrm_nom wrm_crn/4];

% Given Power for 6 operating points
p_max_pk = 100*1000;                     % Watts
p_max_ct = 55*1000;                      % Watts     
p_nom_ct = 27.5*1000;                    % Watts
P = [p_max_pk p_max_ct p_max_ct p_max_ct/3 p_max_ct/3 p_max_pk/4];

% Torque requirement for 6 operating points
Te = P./WRM;

% design requirements------------------------------------------------------
D.vdc  = 1000;                       % dc voltage (V)
D.wrm  = WRM;                        % mechanical rotor speed 
D.Te   = Te;                         % req. torque in motor mode(Nm)
D.rrs  = 2*cm;                       % shaft radius
D.kpf  = 0.5;                        % packing factor                                     
D.leo  = 1*cm;                       % end winding offset
D.vfs  = 2;                          % forward semiconductor drop (V)
D.J    = 20;                         % number of rotor positions
D.mlim  = 40;                        % mass limit (kg)
D.Pllim = 6000;                      % loss limit (W)
D.atar  = 10;                        % maximum tooth aspect ratio
D.aso = 1.5;                         % slot opening factor
D.nspp=2;                            % slots per pole per phase
D.km=0.75;                           % multiple of magnet Hci to get Hlim
D.lfs=0.0;                           % length of front shaft (m) 
D.lbs=0.0;                           % length of back shaft (m) 