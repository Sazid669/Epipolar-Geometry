function [ cam1_fig,cam2_fig] = mvg_show_epipoles(cam1_fig, cam2_fig,ep1,ep2)
%SHOW_EPIPOLAR_GEOM Display epipoles
% Show epipole on image 1
set(0, 'CurrentFigure', cam1_fig)
hold on;
scatter(ep1(1),ep1(2),30,'g','LineWidth',1.5);
% Show lines on image 2
set(0, 'CurrentFigure', cam2_fig)
hold on;
scatter(ep2(1),ep2(2),30,'g','LineWidth',1.5);
end

