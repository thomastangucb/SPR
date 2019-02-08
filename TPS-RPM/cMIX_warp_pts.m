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
% 3 % %%% cMIX_warp_pts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%
function [vx] = cMIX_warp_pts (trans_type, x, z, c_tps, d_tps, w, sigma_kernel);

switch (trans_type)
  case 'tps'
    vx = ctps_warp_pts (x, z, c_tps, d_tps); 
  case 'rbf'
    vx = crbf_warp_pts (x, z, w, sigma_kernel);
  % case 'gtm_tps'
  %  vx = cgtm_warp_pts ('tps_style', x, z, w, 0);
  otherwise; 
    disp ('# ERROR #: cMIX_warp_pts -- wrong input!');
end
