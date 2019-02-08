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


% ---------------------------------------------------------------
% ctps_plot_grid_gen.m
% ---------------------------------------------------------------
% Generate grid points for displaying TPS deformation.
%
% Usage:
% [grid_pts, controls] = ctps_plot_grid_gen (x);
% [grid_pts, controls] = ctps_plot_grid_gen (x, resolution, resolution_grid);
%
%         *** "controls" are for: ctps_plot_gridbox.
% 
% 01/26/00

function [grid_pts, controls] = ctps_plot_grid_gen (x, resolution, resolution_grid);

% check input:
% ------------
if (nargin == 1);      % input (x), set the other 2.
  resolution = 4;
  resolution_grid = 3;
elseif (nargin == 3);  % input (x, resolution, resolution_grid);
  ;
else
  disp ('# ERROR #: ctps_plot_grid_gen -- wrong input!');
  help ctps_plot_grid_gen;
end;



% set grid range:
% ---------------
xrange = [min(x(:,1)), max(x(:,1))];
yrange = [min(x(:,2)), max(x(:,2))];

% expand a little bit:
expand_ratio = 5;
xrange (1) = xrange (1) - (xrange(2)-xrange(1))/expand_ratio;
xrange (2) = xrange (2) + (xrange(2)-xrange(1))/expand_ratio;
yrange (1) = yrange (1) - (yrange(2)-yrange(1))/expand_ratio;
yrange (2) = yrange (2) + (yrange(2)-yrange(1))/expand_ratio;


% Generate the grid points:
% -------------------------
[grid_pts, rows, cols, points_row, points_col] = cgrid_generate ...
    (xrange(1), xrange(2), yrange(1), yrange(2), resolution, resolution_grid);


controls    = zeros (4,1);
controls(1) = rows;
controls(2) = cols;
controls(3) = points_row;
controls(4) = points_col;




%%%%%
% 1 % %%% cgrid_generate.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%
% Generate grid.
%
% Input: x/y       -- one of data sets.
%
% Output: grid_pts  -- grid points, same format as x/y.
%         rows,cols -- grid dimensionality.
%         points_row, points_cols -- points along each
%                                    row and col.
% ------------------------------------------------------------
% Last modified: 09/27/99

function [grid_pts, rows, cols, points_row, points_col] = ...
    cgrid_generate (xrange1, xrange2, yrange1, yrange2, ...
                    resolution, resolution_grid);

xrange = [xrange1, xrange2];
yrange = [yrange1, yrange2];

% a: grid square size.
a = min(xrange(2)-xrange(1), yrange(2)-yrange(1)) / resolution;
grid_step = a / resolution_grid;

rows = ceil((yrange(2)-yrange(1)) / a + 1);
cols = ceil((xrange(2)-xrange(1)) / a + 1);

yrange(2) = yrange(1) + (rows-1)*resolution_grid*grid_step;
xrange(2) = xrange(1) + (cols-1)*resolution_grid*grid_step;

%keyboard

grid_pts = [];
% points_row = floor( (xrange(2)-xrange(1))/grid_step ); 
points_row = (cols-1) * resolution_grid + 1; % two ending points.
for i=1:rows
  tmp_row = [[xrange(1):grid_step:xrange(2)]', ...
	ones(points_row,1) * (i-1) * a + yrange(1)];
  grid_pts = [grid_pts; tmp_row];
end;


% points_col = floor((yrange(2)-yrange(1))/grid_step ); 
points_col = (rows-1) * resolution_grid + 1; % two ending points.
for j=1:cols
  tmp_col = [ones(points_col,1) * (j-1) * a + xrange(1), ...
	[yrange(1):grid_step:yrange(2)]'];
  grid_pts = [grid_pts; tmp_col];
end;

