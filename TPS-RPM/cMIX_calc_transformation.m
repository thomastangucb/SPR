% Robust Point Matching (RPM) Demo (version 20000427):
% ----------------------------------------------------
% Copyright (C) 2000 Haili Chui, Anand Rangarajan
% 
% Authors: Haili Chui and Anand Rangarajan
% Date:    04/27/2000
% 
% Contact Information:
%
% Haili Chui:		chui@noodle.med.yale.edu
% Anand Rangarajan:	anand@noodle.med.yale.edu
% 
% Terms:	  
% 
% The source code (M-files) are provided under the
% terms of the GNU General Public License with an explicit
% clause permitting the execution of the M-files from within
% a MATLAB environment. See the LICENSE file for details.
%
%




%%%%%
% 2 % %%% cMIX_calc_transformation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%
% 
function [c_tps, d_tps, w] = cMIX_calc_transformation (transformation_type, ...
    lamda1, lamda2, sigma_kernel, x, vy, z);

c_tps = [];
d_tps = [];
w     = [];

switch (transformation_type)
  case 'tps'
    [c_tps,d_tps]  = ctps_gen (x, vy, lamda1, lamda2);
    % [c_tps, d_tps] = ctps_generate_cd_regularized (lamda1, lamda2, x, vy);  
  case 'rbf'
    [phi, w] = crbf_gen (x, vy, z, lamda1, lamda2, sigma_kernel);
  %case 'gtm_tps'
  %  [phi, w] = cgtm_calc_w ('tps_style', x, vy, z, lamda1, lamda2, 0);
  otherwise; 
    disp ('# ERROR #: cMIX_calc_transformation -- wrong input!');
end



