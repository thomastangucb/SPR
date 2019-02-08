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
% cplot.m
% ------------------------------------------------------------------- 
% Plot points.
%
% Usage: 
% [] = cplot (x);
% [] = cplot (x,y);
% [] = cplot (x,y,z);
%
% [] = cplot (x, marker_str, marker_size);
% [] = cplot (x, xmarker_str, xmarker_size, y, ymarker_str, y_marker_size);
%
% 02/01/00

function [] = cplot (in1,in2,in3,in4,in5,in6);

% check input:
if (nargin == 1) % -------------------------------- (x) ---
  x = in1; xmarker_str = 'go'; xmarker_size = 6;
  y = [];  ymarker_str = 'r+'; ymarker_size = 6;
  z = [];  zmarker_str = 'bo'; zmarker_size = 6;
  
elseif (nargin == 2) % -------------------------- (x,y) ---
  x = in1; xmarker_str = 'go'; xmarker_size = 6;
  y = in2; ymarker_str = 'r+'; ymarker_size = 6;
  z = [];  zmarker_str = 'bo'; zmarker_size = 6;
  
elseif (nargin == 3) & (~isstr(in2)) % -------- (x,y,z) ---
  x = in1; xmarker_str = 'go'; xmarker_size = 6;
  y = in2; ymarker_str = 'r+'; ymarker_size = 6;
  z = in3; zmarker_str = 'bo'; zmarker_size = 6;
  
elseif (nargin == 3) & (isstr(in2)) % ---- (x, 'go', 3) ---
  x = in1; xmarker_str = in2; xmarker_size = in3;
  y = [];
  z = [];
  
elseif (nargin == 6) % -------- (x, 'go, 3, y, 'r+', 3) ---
  x = in1; xmarker_str = in2; xmarker_size = in3;
  y = in4; ymarker_str = in5; ymarker_size = in6;
  z = [];
  
else
  disp ('# ERROR #: cplot -- wrong input!');
  help cplot; return;
end;

% plot x:
[n, dim] = size(x);
if (n >= 1); cplot_1pointset (x, xmarker_str, xmarker_size); end;
hold on;

% plot y:
[n, dim] = size(y);
if (n >= 1); cplot_1pointset (y, ymarker_str, ymarker_size); end;

% plot z:
[n, dim] = size(z);
if (n >= 1); cplot_1pointset (z, zmarker_str, zmarker_size); end;
hold off;





%%%%%
% 1 % %%% cplot_1pointset %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%
% Plot one point set.
% 
% Usage:
% [] = cplot_1pointset (x, xmarker_str, xmarker_size);
%
% 02/01/00

function [] = cplot_1pointset (x, xmarker_str, xmarker_size);

[n,dim] = size(x);
if (dim == 2)
  h = plot (x(:,1), x(:,2), xmarker_str, 'markersize', xmarker_size); axis('equal');
  hold on;
elseif (dim == 3)
  h = plot3 (x(:,1), x(:,2), x(:,3), xmarker_str, 'markersize', xmarker_size); axis('equal');
  axis('equal'); set (gca, 'box', 'on');
  hold on;
end;
