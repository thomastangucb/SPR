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


% ---------------------------------------------------------
% ctps_plot_grid.m
% ---------------------------------------------------------
% Display TPS deformed grid.
%
% Usage: 
% [] = ctps_plot_grid (x,y,c,d);
% [] = ctps_plot_grid (x,y,lamda);
% [] = ctps_plot_grid (x,y,c,d,  resolution,resolution_grid); 
% [] = ctps_plot_grid (x,y,lamda,resolution,resolution_grid); 
%
% for (x,y,c,d)
% x -- TPS basis points.
% y -- points to be warped.
%
% for (x,y,lamda)
% x -- TPS basis points and points to be warped.
% y -- target points.
%
% 01/26/00


function [] = ctps_plot_grid_simple (x,y,c,d,...
    resolution,resolution_grid);

% t = [7,6; 10,9; 1 1; 1 15; 15 1; 15 15];
% y = [11,7; 7,10; 1 1; 1 15; 15 1; 15 15];
% [c,d] = ctps_gen (t, y, 1);

% check input.
% ------------
if (nargin == 3);     % input (x,y,lamda);
  lamda = c;
  [c,d] = ctps_gen (x,y,lamda);
  tmp = x;  y = x; x = tmp;

  resolution      = 4;
  resolution_grid = 3;

elseif (nargin == 4); % input (x,y,c,d);
  resolution      = 4;
  resolution_grid = 3;

elseif (nargin == 5); % input (x,y,lamda,resolution,resolution_grid); 
  lamda           = c;
  resolution      = d;
  resolution_grid = resolution;
  [c,d]           = ctps_gen (x,y,lamda);
  tmp = x;  y = x; x = tmp;
  
elseif (nargin == 6);  % input (x,y,c,d,resolution,resolution_grid); 
  % do nothing.

else
  disp ('# ERROR #: ctps_plot_grid -- wrong input!');
  help ctps_plot_grid; return;
end;



% generate grid points. 
% ---------------------
[grid_pts, controls] = ctps_plot_grid_gen (x, resolution, resolution_grid);


% Warp the grid:
% --------------
[grid_new] = ctps_warp_pts (grid_pts, x, c, d); 


% Plot:
% -----
ori_color = ones(1,3) * 0.7;
new_color = 'b';

ctps_plot_gridbox (1, grid_pts, controls, ori_color, ':'); hold on;
ctps_plot_gridbox (1, grid_new, controls, 'b','-'); hold on;
%axis('equal'); axis ('off');


