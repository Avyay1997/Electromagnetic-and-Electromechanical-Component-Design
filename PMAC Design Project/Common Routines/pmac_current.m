function [I] = pmac_current(I)
% pmac_current calculates the rms current and phase based on q- and d-
%              axis current
%
% Call:
% I=pmac_current(I)
%
% Inputs:
% I        = structure of current related variables. input fields are
%  I.iqsr  = q-axis current vector - one element per operating point (A)
%  I.idsr  = d-axis current vector - one elemnt per operating point (A)
% 
% I        = structure of current related variables. output fields are
%  I.Is    = vector of RMS phase current amplitude - one element per 
%            operating point (A)
%  I.phii  = vector of current phase angles - one element per operating 
%            point (rad)
%
% Written by:
% S.D. Sudhoff
% School of Electrical and Computer Engineering
% 1285 Electrical Engineering Building
% West Lafayette, IN 47907-1285
% Phone: 765-494-3246
% Fax: 765-494-0676
% E-mail: sudhoff@ecn.purdue.edu
               
I.Is=sqrt((I.iqsr.^2+I.idsr.^2)/2);
I.phii=angle(I.iqsr - 1j*I.idsr);