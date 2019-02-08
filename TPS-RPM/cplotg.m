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


% ----------------------------------------------------------------
% cplot_2g.m
% ----------------------------------------------------------------
% Plot links between 2 corresponding point sets.
%
% Usage: 
% [] = cplotg (x, y);                % m = eye.
% [] = cplotg (x, y, m);             % thr = 0.
% [] = cplotg (x, y, m, threshold);
% [] = cplotg (x, y, m, threshold, color_str);
% [] = cplotg (x, y, m, threshold, color_str, marker_str);
% 
% 02/01/00

function [] = cplotg(x, y, m, threshold, color_str, marker_str);

% check input:
% ------------
[n, dim] = size (x);

if (nargin == 2)
  m         = eye (n,n);
  threshold = 0;
  gcolor_str = 'y:';

elseif (nargin == 3)
  threshold = 0;
  gcolor_str = 'y:';

elseif (nargin == 4)
  gcolor_str = 'y:';

elseif (nargin == 5)
  gcolor_str = color_str;
  
elseif (nargin == 6)
  ;
else
  disp ('# ERROR #: cplotg -- wrong input!');
  help cplotg; return;
end;

%keyboard

switch (dim)
  
  case 2 % --------------------------------------- 2D ---
    
    % normal plot:
    if (nargin < 6)
      % Reformat data:
      % --------------
      xy = [x; y];
      
      [siz1, temp] = size(x);
      [siz2, temp] = size(y);
      
      [indexi, indexj] = find ( m > threshold );
      index = indexi + (indexj-1) * siz1;
      
      msp = zeros(siz1,siz2);
      msp (index) = 1;
      
      madj = [zeros(siz1), msp; 
	msp',zeros(siz2)];
      
      
      % Plot:
      % -----
      axis ('equal'); axis('off'); 
      gplot (madj, xy, gcolor_str);
  
    % otherwise, i want to change the color:
    else 
      [siz1, temp] = size(x);
      [siz2, temp] = size(y);
      
      [indexi, indexj] = find ( m > threshold );
      n = length(indexi);    
      
      % Plot:
      % -----
      for i=1:n
	hold on; 
	tmp = [x(indexi(i),:); y(indexj(i),:)];
	plot (tmp(:,1), tmp(:,2), ... 
	    'color', color_str, 'linestyle', marker_str);
	axis ('equal'); 
      end;
    end;
    
  case 3 % ----------------------------------------- 3D ---
    [siz1, temp] = size(x);
    [siz2, temp] = size(y);
    
    [indexi, indexj] = find ( m > threshold );
    n = length(indexi);    
    
    % Plot:
    % -----
    for i=1:n
      hold on;
      tmp = [x(indexi(i),:); y(indexj(i),:)];
      plot3 (tmp(:,1), tmp(:,2), tmp(:,3), ... 
	  'color', color_str, 'linestyle', marker_str);
    end;
  otherwise;
end;


