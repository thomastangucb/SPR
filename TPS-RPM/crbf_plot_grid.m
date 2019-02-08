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


% ------------------------------------------------------------------- 
% crbf_plot_grid.m
% ------------------------------------------------------------------- 
% Plot RBF warpped grid.
%
% Usage: 
% [] = crbf_plot_grid (x,z,w,sigma_kernel);
%      x is the point set to be warped.
%
% 01/26/00

function [] = crbf_plot_grid (x,z,w,sigma_kernel);

if (nargin ~= 4)
  disp ('# ERROR #: crbf_plot_grid -- wrong input!');
  help crbf_plot_grid; return;
end;


% generate grid:
src = x;

[grid_pts, controls] = ctps_plot_grid_gen (src);

% warp grid:
grid_pts1 = crbf_warp_pts (grid_pts,z,w,sigma_kernel);

% display:
ori_color = ones(1,3) * 0.7;

ctps_plot_gridbox (1, grid_pts,  controls, ori_color, ':');
ctps_plot_gridbox (1, grid_pts1, controls, 'b','-');
axis('equal'); axis('off');
