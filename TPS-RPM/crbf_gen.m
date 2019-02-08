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
% crbf_gen.m
% ------------------------------------------------------------------- 
% Generate Radial-Basis-Function spline parameters.
%
% Usage: 
% Usage: 
% [phi]   = crbf_gen (x, z, sigma_kernel);
% [phi,w] = crbf_gen (x, y, z, lamda1, sigma_kernel);
% [phi,w] = crbf_gen (x, y, z, lamda1, lamda2, sigma_kernel);
% 
% x -- pts to be warpped.
% y -- target pts.
% z -- basis pts.
%
% 01/26/00

function [o1,o2] = crbf_gen (x,y,z, lamda1, lamda2, sigma_kernel);


% check input:
% ------------

% --- [phi]   = gen_crbf (x, z, sigma_kernel) ----------------------------
if (nargin == 3); 
  sigma_kernel = z;
  z            = y;

  [L,dim] = size(x);
  P       = length(z);

  [phi] = crbf_kernel (x, z, sigma_kernel);
  
  o1 = phi;

% --- [phi,w] = crbf_gen (x, y, z, lamda1, sigma_kernel) -----------------
elseif (nargin == 5);
  sigma_kernel = lamda2;
  lamda2       = 0;

  [L,dim] = size(x);
  P       = length(z);

  [phi] = crbf_kernel (x, z, sigma_kernel);

  G       = eye(L,L);                   % assume no outlier in "x".
  
  Iw                                 = zeros (P+dim+1, P+dim+1); 
  Iw (1:P,1:P)                       = lamda1 * eye(P,P);
  Iaffine                            = zeros (P+dim+1, P+dim+1); 
  Iaffine (P+1:P+dim+1, P+1:P+dim+1) = lamda2 * eye (dim+1,dim+1); 
  
  w = inv(phi' * G * phi + (Iw + Iaffine)) * phi' * y;

  o1 = phi;
  o2 = w;
  
  
  
% --- [phi,w] = crbf_gen (x, y, z, lamda1, lamda2, sigma_kernel) ---------
elseif (nargin == 6);
  [L,dim] = size(x);
  P       = length(z);

  [phi] = crbf_kernel (x, z, sigma_kernel);

  G       = eye(L,L);                   % assume no outlier in "x".
  
  Iw                                 = zeros (P+dim+1, P+dim+1); 
  Iw (1:P,1:P)                       = lamda1 * eye(P,P);
  Iaffine                            = zeros (P+dim+1, P+dim+1); 
  Iaffine (P+1:P+dim+1, P+1:P+dim+1) = lamda2 * eye (dim+1,dim+1); 
  
  w = inv(phi' * G * phi + (Iw + Iaffine)) * phi' * y;

  o1 = phi;
  o2 = w;
  
  
else
  disp ('# ERROR #: crbf_gen -- wrong input!');
  help crbf_gen;
end;










%%%%%
% 1 % %%% crbf_kernel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%
%
% Calc. BRF kernel "phi".
%
% Usage:
% [phi] = cgtm_calc_kernel (x, z, sigma_kernel)
%
% 01/26/00

function [phi] = crbf_kernel (x, z, sigma_kernel)

[L,dim] = size(x);
P       = length(z);

phi = zeros(L,P);

for it_dim=1:dim
  tmp = x(:,it_dim) * ones(1,P) - ones(L,1) * z(:,it_dim)';
  tmp = tmp .* tmp;
  phi = phi + tmp;
end;

if dim == 2
  phi = exp(-phi/(sigma_kernel^2)); % gaussian basis fn.
  
else
  phi = - sqrt(phi);
end;

phi = [phi, x, ones(L,1)];   % x, ones: affine part.





