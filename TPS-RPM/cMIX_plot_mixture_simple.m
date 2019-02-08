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


% cMIX_plot_mixture_simple.m
% ------------------------------------------------------------------- 
% Plot the clusters as bunch of ecllipses.
%
% Usage: [] = cMIX_plot_mixture_simple (vx, T);
% ------------------------------------------------------------------- 
% Last modified: 02/01/00

function [] = cMIX_plot_mixture_simple (vx, T, color_str);

% check input:
% ------------
if (nargin == 2)
  color_str = 'b-'; % default.
end;

[c,dim] = size (vx);

a = sqrt(T);
b = sqrt(T);

% Generate the ecllips:
% ---------------------
step_theta = 20;
tmp_theta  = [1:step_theta:360+step_theta]'/360*2*pi;
tmp_pts    = [a.*cos(tmp_theta), b.*sin(tmp_theta)];
n          = length(tmp_theta);

for i=1:c
  if (dim == 2)
    % Draw the ecllips:
    % -----------------
    plot (tmp_pts(:,1) + vx(i,1), tmp_pts(:,2)+vx(i,2), color_str); hold on;
  else
    % disp ('plot3');
    plot3 (tmp_pts(:,1)+vx(i,1), tmp_pts(:,2)+vx(i,2), zeros(n,1)+vx(i,3), color_str); hold on;
    %plot3 (zeros(n,1)+vx(i,1), tmp_pts(:,1)+vx(i,2), tmp_pts(:,2)+vx(i,3), color_str); hold on;
    %plot3 (tmp_pts(:,1)+vx(i,1), zeros(n,1)+vx(i,2), tmp_pts(:,2)+vx(i,3), color_str); hold on;

    %[X,Y,Z] = sphere(8);
    
    %X1 = X*T + vx(i,1);
    %Y1 = Y*T + vx(i,2);
    %Z1 = Z*T + vx(i,3); hold on;
    %mesh (X1, Y1, Z1); caxis([19 19.000000001]);hidden off;
  end;
end;

