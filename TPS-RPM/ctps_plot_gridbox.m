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
% ctps_plot_gridbox.m
% ------------------------------------------------------------------- 
% Plot grid box for TPS demonstration.
%
% Usage: 
% [] = cplot_gridbox (grid_method, grid, controls);
% 
% grid_method: 0 -- pts
%              1 -- pts + linking lines.
% controls: (rows, cols, points_row, points_col).
%
% 01/26/00

function [] = ctps_plot_gridbox (grid_method, grid, controls, ...
    marker_color, marker_type);

rows       = controls(1);
cols       = controls(2);
points_row = controls(3);
points_col = controls(4);



switch (grid_method)
  
  case 0
    plot (grid(:,1), grid(:,2), '.', 'color', marker_color,'markersize', 5);
    
  case 1
    for i=1:rows
      tmp = grid ( (i-1)*points_row+1:i*points_row,:);
      plot (tmp(:,1), tmp(:,2), 'color', marker_color, ...
	  'linestyle', marker_type);%, ...
	  %'erasemode', 'background');
      hold on;
    end;
    
    start_index = rows * points_row;
    for j=1:cols
      tmp = grid ( (j-1)*points_col + start_index + 1 : j*points_col + ...
	  start_index, :);
      plot (tmp(:,1), tmp(:,2), 'color', marker_color, ...
	  'linestyle', marker_type);%, ...
%	  'erasemode', 'background');
      hold on;
    end;
  otherwise;
    disp ('ERROR: cplot_gridbox -- wrong input parameters');
end;

