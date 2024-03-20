function [ lm1,lm2,l_coef_1,l_coef_2 ] = mvg_compute_epipolar_geom_modif( cam1_p2d,cam2_p2d,F)
%SHOW_EPIPOLAR_GEOM Computing epipolar line cooeficients 
if size(cam1_p2d,1) == 2
    cam1_p2d = [cam1_p2d; ones(1,size(cam1_p2d,2))];
end
if size(cam2_p2d,1) == 2
    cam2_p2d = [cam2_p2d; ones(1,size(cam2_p2d,2))];
end

% Epiploar line on image plane 2
lm2 = F * cam1_p2d;
% Epiploar line on image plane 1
lm1 = F' * cam2_p2d;
%Line coeff on image plane 1
l_coef_1(1,:) = -lm1(1,:)./lm1(2,:);
l_coef_1(2,:) = -lm1(3,:)./lm1(2,:);
%Line coeff on image plane 2
l_coef_2(1,:) = -lm2(1,:)./lm2(2,:);
l_coef_2(2,:) = -lm2(3,:)./lm2(2,:);
end

