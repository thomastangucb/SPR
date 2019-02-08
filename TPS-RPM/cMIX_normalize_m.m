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


% cMIX_normalize_m.m
% ------------------------------------------------------------------- 
% Double normalization of m (with outlier col/row).
%
% Input:  m, m_outlier_col, m_outlier_row -- m entries.
% Output: m, m_outlier_col, m_outlier_row -- normalized.
% ------------------------------------------------------------------- 
% Last modified: 10/11/99

function [m, m_outlier_col, m_outlier_row] = cMIX_normalize_m (m, ...
    m_outlier_col, m_outlier_row);

% Parameters:
norm_threshold = 0.05;
norm_maxit     = 10;

[xmax, ymax] = size(m);
  
norm_it = 0;
while (1 > 0)

  % --- Row normalization --------------------------------------------
  sx = sum(m')' + m_outlier_col;
  m  = m ./ (sx * ones(1,ymax));
  m_outlier_col = m_outlier_col ./sx;
  
  % --- Column normalization -----------------------------------------
  sy = sum(m) + m_outlier_row;
  m  = m ./ (ones(xmax,1)*sy);
  m_outlier_row = m_outlier_row ./sy;
  
  % time to quit?
  err = ((sx-1)'*(sx-1) + (sy-1)*(sy-1)')/(length(sx)+length(sy));
  if err < (norm_threshold .* norm_threshold); break, end
  % run out of time:
  norm_it = norm_it + 1;
  if norm_it >= norm_maxit; break; end
  
end
