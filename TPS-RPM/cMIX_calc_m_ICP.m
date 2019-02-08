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


% ------------------------------------------------------------
% cMIX_calc_m_ICP.m:
% ------------------------------------------------------------
% Calc. two-way  ICP m.
% 
% Usage:
% [m, dist_threshold] = cMIX_calc_m_ICP (vx,y);
% [m, dist_threshold] = cMIX_calc_m_ICP (vx,y,k_sigma);
% [m, dist_threshold] = cMIX_calc_m_ICP (vx,y,k_sigma,which_way);
%
% Notes: which_way -- specify one way ICP m.
%                0 -- both way (default).
%                1 -- x to y.
%                2 -- y to x.
% 02/03/00

function [m, dist_threshold] = cMIX_calc_m_ICP (vx,y,k_sigma,which_way);

% check input:
if (nargin == 2)
  k_sigma = 0;   % default, no outlier rejection.
  which_way = 0;

elseif (nargin == 3)
  which_way = 0;

elseif (nargin == 4)
  
else
  disp ('# ERROR #: cMIX_calc_m_ICP -- wrong input !');
  help cMIX_calc_m_ICP; return;
end;

  
if k_sigma == 0; 
  dist_threshold_flag = 0;
  dist_threshold      = 1e10; 
else
  dist_threshold_flag = 1;
  dist_threshold      = 0;
end; % no outlier.

[siz1, dim]  = size(vx); xmax = siz1;
[siz2, temp] = size(y); ymax = siz2;

% Find nearest neighbour for each pt in x:
% ----------------------------------------
[M1, dist_x] = cICP_findneighbours (vx, y);
[M2, dist_y] = cICP_findneighbours (y, vx);

if dist_threshold_flag ~= 0
  dist    = [dist_x; dist_y];
  n       = length (dist);
  mean_x  = sum (dist) / n;
  sx      = std(dist); % sqrt ( 1 / (n-1) * sum ( (dist - mean_x) .* (dist - mean_x)));
  
  dist_threshold = mean_x + k_sigma * sx;
end;


% this update thing of threshold doesn't work !!!
% -----------------------------------------------
% xdist_sum = sum (xdistances);
% dist_threshold = xdist_sum / (siz1/frac) * 3;
m1 = zeros(xmax,ymax); m2 = m1; 
for i=1:xmax
  if dist_x(i) > dist_threshold; M1 (i) = -1; else; m1(i, M1(i)) = 1; end;
end;
for j=1:ymax
  if dist_y(j) > dist_threshold; M2 (j) = -1; else; m2(M2(j), j) = 1; end;
end;

if (which_way == 0)
  m = (m1 +m2)/2;
elseif (which_way == 1)
  m = m1;
elseif (which_way == 2)
  m = m2;
end;







%%%%%
% 1 % %%% cICP_fnneighbours %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%
% 
% Find nearest neighbour in template t for each point in x. 3d/2d
%
% Usage: [M, distances] = cICP_findneighbours (x, t);
% Notes: 
%         M -- M(1) = 10, y(10) is nearest from x(1);
%         distance -- distance (x(i)-t(j))^2
% 
% 02/01/00

function [M, distances] = cICP_fnneighbours (x, t);

[m, dim] = size(x);
[n, dim] = size(t);
M = zeros (m,1);

% |x-t| matrices:
% ---------------
xttmp = zeros (n, m);
for i=1:dim
  xtmp = ones(n,1) * x(:,i)';
  ttmp = t(:,i)  * ones(1,m); 
  xttmp = xttmp + (xtmp - ttmp) .* (xtmp - ttmp);
end;

% M + min_dist list:
% ------------------
[min_dist, min_index] = min(xttmp);
distances = (sqrt(min_dist))';
M         = min_index';




