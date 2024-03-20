function [ cam1_fig,cam2_fig] = mvg_show_epipolar_lines(cam1_fig, cam2_fig,l_coef_1,l_coef_2,margins,clr)
%SHOW_EPIPOLAR_GEOM Display epipolar lines

%Image plane 1
%Get coordinates of start and end point of line 
l_x_1 = zeros(2,size(l_coef_1,2));
l_x_1(1,:) = margins(1,1);     % Start x of line
l_x_1(2,:) = margins(2,1);   % End x of line 
l_y_1(1,:) = l_x_1(1,:).*l_coef_1(1,:)+l_coef_1(2,:);    % Start y of line
l_y_1(2,:) = l_x_1(2,:).*l_coef_1(1,:)+l_coef_1(2,:);    % End y of line

%Image plane 2
%Get coordinates of start and end point of line 
l_x_2 = zeros(2,size(l_coef_2,2));
l_x_2(1,:) = margins(1,2);     % Start x of line
l_x_2(2,:) = margins(2,2);   % End x of line 
l_y_2(1,:) = l_x_2(1,:).*l_coef_2(1,:)+l_coef_2(2,:);    % Start y of line
l_y_2(2,:) = l_x_2(2,:).*l_coef_2(1,:)+l_coef_2(2,:);    % End y of line

% Show lines and epipole on image 1
set(0, 'CurrentFigure', cam1_fig)
hold on;
line(l_x_1,l_y_1,'Color',clr);
axis auto
axis equal
% Show lines on image 2
set(0, 'CurrentFigure', cam2_fig)
hold on;
line(l_x_2,l_y_2,'Color',clr);
axis auto
axis equal
end

