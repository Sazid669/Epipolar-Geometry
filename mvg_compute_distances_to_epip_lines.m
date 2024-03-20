function [ absDistVec1, absDistVec2] = mvg_compute_distances_to_epip_lines( cam1_p2d,cam2_p2d,F)
% Computes two vectors of distances of points to corresponding epipolar lines 

[ lm1,lm2,l_coef_1,l_coef_2 ] = mvg_compute_epipolar_geom_modif( cam1_p2d,cam2_p2d,F);

[ res_d1, absDistVec1  ] = mvg_compute_distance_point_line(lm1,cam1_p2d);
[ res_d1, absDistVec2  ] = mvg_compute_distance_point_line(lm2,cam2_p2d);


return;

