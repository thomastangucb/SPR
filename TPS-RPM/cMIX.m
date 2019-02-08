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


% ----------------------------------------------------------------------
% cMIX.m
% ---------------------------------------------------------------------- 
% Gaussian Mixture Robust Point Matching.
% 
% Usage: 
%
% [c,d,m] = cMIX (x, y, frac, Tinit, Tfinalfac);
%
% Optional ones:
%
% [c,d] = cMIX (x, y, frac, Tinit, Tfinalfac);
% [c,d] = cMIX (x, y, frac, Tinit, Tfinalfac, disp_flag);
% [c,d] = cMIX (x, y, frac, Tinit, Tfinalfac, disp_flag, m_method);
%
% [w]   = cMIX (x, y, z, sigma, frac, Tinit, Tfinalfac);
% [w]   = cMIX (x, y, z, sigma, frac, Tinit, Tfinalfac, disp_flag);
% [w]   = cMIX (x, y, z, sigma, frac, Tinit, Tfinalfac, disp_flag, m_method, 
% icp_sigma);
%
% Notes: the program will set transformation type automatically.
%        ['tps', 'rbf'].
%
%        'icp' -- 'icp0', 'icp3', 'icp5'. --> to set k_sigma.
% 
% 04/27/00

function [o1,o2,o3 , o4] =  cMIX (in1,in2,in3,in4,in5,in6,in7,in8,in9);

%fig(1); clf; whitebg('k'); set(gcf,'color',[0 0 0]);
% set(gcf,'DoubleBuffer','on')

% Init control parameters:
% ------------------------
perT_maxit  = 5;
relax_maxit = 1;
anneal_rate = 0.93;

lamda1_init = 1;
lamda2_init = 0.01;

% default:
% --------
disp_flag = 1;
m_method  = 'mixture';

debug_flag = 0;


% check input:
% ------------

% --- [c,d] = cMIX (x, y, frac, Tinit, Tfinalfac) ---------------------------
if (nargin == 5) 
  x          = in1;
  y          = in2;
  frac       = in3;
  T_init     = in4;
  T_finalfac = in5;

  trans_type = 'tps';  
  
  z         = x;
  sigma     = 1;
  
% --- [c,d] = cMIX (x, y, frac, Tinit, Tfinalfac, disp_flag) --------------
elseif (nargin == 6) & (length(in3) == 1)
  x          = in1;
  y          = in2;
  frac       = in3;
  T_init     = in4;
  T_finalfac = in5;
  disp_flag  = in6;

  trans_type = 'tps';  

  z          = x;
  sigma      = 1;

% --- [c,d] = cMIX (x, y, frac, Tinit, Tfinalfac, disp_flag, m_method) ------
elseif (nargin == 7) & (length(in3) == 1)
  x          = in1;
  y          = in2;
  frac       = in3;
  T_init     = in4;
  T_finalfac = in5;
  disp_flag  = in6;
  m_method   = in7;
  
  trans_type = 'tps';  

  z          = x;
  sigma      = 1;


% --- [w]   = cMIX (x, y, z, sigma, frac, Tinit, Tfinalfac) -----------------
elseif (nargin == 7) & (length(in3) > 1)
  x          = in1;
  y          = in2;
  z          = in3;
  sigma      = in4;
  frac       = in5;
  T_init     = in6;
  T_finalfac = in7;

  trans_type = 'rbf';  

  disp_flag = 1;
  
  
% --- [w]   = cMIX (x, y, z, sigma, frac, Tinit, Tfinalfac, disp_flag) -----
elseif (nargin == 8) & (length(in3) > 1)
  x          = in1;
  y          = in2;
  z          = in3;
  sigma      = in4;
  frac       = in5;
  T_init     = in6;
  T_finalfac = in7;
  disp_flag  = in8;
  
  trans_type = 'rbf';  

% --- [w]   = cMIX (x,y,z,sigma,frac,Tinit,Tfinalfac,disp_flag,m_method); ---
elseif (nargin == 9) & (length(in3) > 1)
  x          = in1;
  y          = in2;
  z          = in3;
  sigma      = in4;
  frac       = in5;
  T_init     = in6;
  T_finalfac = in7;
  disp_flag  = in8;
  m_method   = in9;
  
  trans_type = 'rbf';  
  
% --- [theta,tx,ty] = cMIX (x, y, frac, Tinit, Tfinalfac) ----------------------
elseif (nargin == 5) & (nargout == 3)
  x          = in1;
  y          = in2;
  frac       = in3;
  T_init     = in4;
  T_finalfac = in5;

  trans_type = 'r+t';  
  theta = 0; t = zeros (2,1); s = 1;

  disp_flag = 1;
  z         = x;
  sigma     = 1;
  
% --- [theta,tx,ty] = cMIX (x, y, frac, Tinit, Tfinalfac, disp_flag) ---------
elseif (nargin == 6) & (nargout == 3)
  x          = in1;
  y          = in2;
  frac       = in3;
  T_init     = in4;
  T_finalfac = in5;
  disp_flag  = in6;

  trans_type = 'r+t';  
  theta = 0; t = zeros (2,1); s = 1;

  z          = x;
  sigma      = 1;
  
% --- [theta,tx,ty] = cMIX (x, y, frac, Tinit, Tfinalfac,disp_flag, m_method); ---
elseif (nargin == 7) & (nargout == 3)
  x          = in1;
  y          = in2;
  frac       = in3;
  T_init     = in4;
  T_finalfac = in5;
  disp_flag  = in6;
  m_method   = in7;

  trans_type = 'r+t';  
  theta = 0; t = zeros (2,1); s = 1;

  z          = x;
  sigma      = 1;
  
else
  disp ('# ERROR #: cMIX -- wrong input!');
  help cMIX; return;
end;

% take care of 'icp' k_sigma stuff.
if (strcmp(m_method(1:3), 'icp'))
  if length(m_method) == 3
    k_sigma  = 0;
    m_method = 'icp';
  else
    k_sigma  = str2num(m_method(4));
    m_method = 'icp';
  end;
else
  k_sigma = 0;
end;


% Init:
% -----

% init x,y,z:
[xmax, dim] = size(x); x = x (1:frac:xmax, :); [xmax, dim] = size(x); 
[ymax, tmp] = size(y); y = y (1:frac:ymax, :); [ymax, tmp] = size(y);  
[zmax, tmp] = size(z); 

if strcmp(trans_type, 'tps')
  z = x;
end;

% init m:
m              = ones (xmax, ymax) ./ (xmax * ymax); 

T0             = max(x(:,1))^2;
moutlier       = 1/sqrt(T0)*exp(-1);       % /xmax *0.001;
% moutlier       = 1/xmax*0.01;
m_outliers_row = ones (1,ymax) * moutlier; 
m_outliers_col = ones (xmax,1) * moutlier; 

% init transformation parameters:
theta = 0; t = zeros (2,1); s = 1;

c_tps = zeros (xmax,dim+1); 
d_tps = eye   (dim+1, dim+1); 

w     = zeros (xmax+dim+1, dim);

% icp: perT_maxit = 1:
if strcmp (m_method(1:3),'icp')
  perT_maxit = 1;
end;

% -------------------------------------------------------------------
% Annealing procedure:
% -------------------------------------------------------------------

T       = T_init;
T_final = T_init / T_finalfac;

vx = x;
vy = y;

it_total = 1;
flag_stop = 0;
while (flag_stop ~= 1)
  
  for i=1:perT_maxit     % repeat at each termperature.

    % Given vx, y, Update m:
    if debug_flag; disp ('calc m ...'); end;
    m = cMIX_calc_m (vx, y, T, m_method, m_outliers_row, m_outliers_col,it_total,k_sigma);

    % Given m, update transformation:
    vy = m * y ./ ( (sum(m'))' * ones(1,dim));
    
    lamda1 = lamda1_init*length(x)*T; 
    lamda2 = lamda2_init*length(x)*T;
    %lamda1 = lamda1_init;
    %lamda2 = lamda2_init;
    
    if debug_flag; disp ('calc c,d ...'); end;
    [c_tps, d_tps, w] = cMIX_calc_transformation (trans_type, lamda1, lamda2, sigma, x, vy, z);

    % w(1:length(z),:) = w(1:length(z),:)*0
    % d_tps = d_tps * 0 + eye(dim+1,dim+1); 
    % c_tps = c_tps *0;
    
    if debug_flag; disp ('calc new x ...'); end;
    [vx] = cMIX_warp_pts (trans_type, x, z, c_tps, d_tps, w, sigma);

  end  % end of iteration/perT
  
  
  T = T * anneal_rate;
  
  % Determine if it's time to stop:
  % -------------------------------
  %  when T <= cluster variance, stop.

  if T < T_final; flag_stop = 1; end;
  


  % Display:
  % -------- 
  it_total = it_total + 1;
  if (disp_flag == 1) |  ((disp_flag == 3) &&(flag_stop == 1)) | ...
	(disp_flag ==2 & mod(it_total,10)==0) | ...
	(it_total == 1)
    figure(1);  clf;
    cMIX_plot_simple (2, x, y, z, vx, m, 1/xmax, T, trans_type, c_tps, d_tps, w, sigma, m_method);  
    str = sprintf ('%s, T = %.4f:\t lamda1: %.4f lamda2: %.4f', m_method,T, lamda1, lamda2);
    disp(str);
    pause(0.1);
  end;

end % end of annealing.



% return outputs:
% ---------------
if strcmp (trans_type, 'tps')
  o1 = c_tps;
  o2 = d_tps;
  o3 = [];
  o4 = vx;
elseif strcmp (trans_type, 'rbf')
  o1 = w;
  o2 = [];
  o3 = [];
elseif strcmp (trans_type, 'r+t')
  o1 = theta;
  o2 = tx;
  o3 = ty;
end;
o3 = m;




%%%%%
% 1 % %%% cMIX_calc_m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%
%
% Update m (correspondence).
%
% Usage:
% [m] = cMIX_calc_m (vx, y, T, 'icp');
% [m] = cMIX_calc_m (vx, y, T, 'mixture');
% [m] = cMIX_calc_m (vx, y, T, 'rpm');
%
% Notes: for "icp", set k_sigma = 0 -- no outlier.
%
% 01/31/00

function [m, m_outliers_row, m_outliers_col] = cMIX_calc_m ...
    (vx, y, T, m_method, m_outliers_row, m_outliers_col, it_total, icp_sigma);

[xmax,dim] = size(vx);
[ymax,dim] = size(y);


% ------------------------------------------------------------------ ICP ---
if strcmp (m_method, 'icp')

  k_sigma = icp_sigma;
  [m, dist_threshold] = cMIX_calc_m_ICP (vx, y, k_sigma);
  m = m + randn(xmax, ymax) * (1/xmax) * 0.001; 
  

% ------------------------------------------------------ one way mixture ---  
elseif strcmp (m_method, 'mixture')
	
  % Given v=tranformed(x), update m:
  y_tmp = zeros (xmax, ymax);
  for it_dim=1:dim
    y_tmp = y_tmp + (vx(:,it_dim) * ones(1,ymax) - ones(xmax,1) * y(:,it_dim)').^2;
  end;  
  
  m_tmp = 1/sqrt(T) .* exp (-y_tmp/T); 
  m_tmp = m_tmp + randn(xmax, ymax) * (1/xmax) * 0.001; 
  
  m = m_tmp;
  
  % normalize accross the outliers as well:
  sy         = sum (m) + m_outliers_row;
  m          = m ./ (ones(xmax,1) * sy); 

  % sx = sum(m')' + m_outliers_col;
  % m2 = m ./ (sx * ones(1,ymax));
  % m = (m+m2)/2;

  
% -------------------------------------------------------- mixture - RPM ---
elseif strcmp (m_method, 'mix-rpm')
	
  % Given v=tranformed(x), update m:
  y_tmp = zeros (xmax, ymax);
  for it_dim=1:dim
    y_tmp = y_tmp + (vx(:,it_dim) * ones(1,ymax) - ones(xmax,1) * y(:,it_dim)').^2;
  end;  
  
  m_tmp = 1/sqrt(T) .* exp (-y_tmp/T); 
  m_tmp = m_tmp + randn(xmax, ymax) * (1/xmax) * 0.001; 
  
  m = m_tmp;

  [m, junk1, junk2] = cMIX_normalize_m (m_tmp, m_outliers_col, m_outliers_row);

  % normalize accross the outliers as well:
  %sy         = sum (m) + m_outliers_row;
  %m          = m ./ (ones(xmax,1) * sy); 

  %sx = sum(m')' + m_outliers_col;
  %m2 = m ./ (sx * ones(1,ymax));
  %m = (m+m2)/2;

  
% --------------------------------------------- RPM, double normalization ---
elseif strcmp (m_method, 'rpm')
  % Given v=tranformed(x), update m:
  y_tmp = zeros (xmax, ymax);
  for it_dim=1:dim
    y_tmp = y_tmp + (vx(:,it_dim) * ones(1,ymax) - ones(xmax,1) * y(:,it_dim)').^2;
  end;  
  
  m_tmp = exp (-y_tmp/T); 
  m_tmp = m_tmp + randn(xmax, ymax) * (1/xmax) * 0.001; 
  
  % double normalization, but keep outlier entries constant.
  moutlier       = 1/xmax * 0.1;
  m_outliers_row = ones (1,ymax) * moutlier; 
  m_outliers_col = ones (xmax,1) * moutlier; 
  
  [m, junk1, junk2] = cMIX_normalize_m (m_tmp, m_outliers_col, m_outliers_row);
  
  
% --------------------------------------------- RPM, double normalization ---
elseif strcmp (m_method, 'rpm-old')
  % Given v=tranformed(x), update m:
  y_tmp = zeros (xmax, ymax);
  for it_dim=1:dim
    y_tmp = y_tmp + (vx(:,it_dim) * ones(1,ymax) - ones(xmax,1) * y(:,it_dim)').^2;
  end;  
  
  m_tmp = exp (-y_tmp/T); 
  m_tmp = m_tmp + randn(xmax, ymax) * (1/xmax) * 0.001; 
  
  % double normalization, also update outlier entries.
  if (it_total == 1)
    moutlier       = 1/xmax * 0.1;
    m_outliers_row = ones (1,ymax) * moutlier; 
    m_outliers_col = ones (xmax,1) * moutlier; 
  end;
  [m, m_outliers_row, m_outliers_col] = cMIX_normalize_m (m_tmp, m_outliers_col, m_outliers_row);
  
else	
  disp ('# ERROR #: cMIX_calc_m -- wrong input!');
end
    





