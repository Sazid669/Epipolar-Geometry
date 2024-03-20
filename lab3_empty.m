% This a empty script to help you move faster on you lab work.
clear all;
close all;
clc;


%% Step 1
% Camera 1
au1 = 100; av1 = 120; uo1 = 128; vo1 = 128;
imageSize = [256 256];
R_c1 = eye(3); % 3x3 Identity matrix
t_c1 = zeros(3,1); % Translation vector [0; 0; 0]

%% Step 2
% Camera 2
au2 = 90; av2 = 110; uo2 = 128; vo2 = 128; 
ax = 0.1; by = pi/4; cz = 0.2; % XYZ EULER 
tx = -1000; ty = 190; tz = 230; 
t_c2 = [tx; ty; tz]; % Translation vector for camera 2
% Construct rotation matrix for camera 2 using XYZ Euler angles
Rx = [1, 0, 0; 0, cos(ax), -sin(ax); 0, sin(ax), cos(ax)]; % Rotation around x-axis
Ry = [cos(by), 0, sin(by); 0, 1, 0; -sin(by), 0, cos(by)]; % Rotation around y-axis
Rz = [cos(cz), -sin(cz), 0; sin(cz), cos(cz), 0; 0, 0, 1]; % Rotation around z-axis
% The rotation matrix for camera 2 is the product of the rotations
R_c2 = Rx * Ry * Rz;


%% STEP 3
% Compute intrinsic matrices and projection matrices

K1 = [au1 0 uo1; 0 av1 vo1; 0 0 1]; % Intrisics matrix for camera 1
wR1c = eye(3);   % rotation of camera 1, from the camera to the world coordinate frame
wt1c = [0 0 0]'; % translation of camera 1, from the camera to the world coordinate frame

% Note: ************** You have to add your own code from here onward ************
K2 = [au2, 0, uo2; 0, av2, vo2; 0, 0, 1]
T_c1w = [R_c1, t_c1; 0, 0, 0, 1]
% Transformation matrix from camera 2 to camera 1
T_c2c1 = [R_c2, t_c2; 0, 0, 0, 1];

% Transformation matrix from camera 2 to world
T_c2w = T_c1w * T_c2c1;

% The inverse of transformation from camera to world is the transformation from world to camera
T_w2c = inv(T_c2w);

% Projection matrix P for camera 1 (since camera 1 is aligned with the world frame)
P1 = K1 * [eye(3), zeros(3,1)];

% Projection matrix P for camera 2
P2 = K2 * T_w2c(1:3,:);

%% STEP 4


% First, calculate the skew-symmetric matrix of t_c2
t_x = [0, -tz, ty; tz, 0, -tx; -ty, tx, 0];

% Then calculate the Fundamental matrix F1
F1 = inv(K2).' * R_c2.' * t_x * inv(K1);
F1 = F1./F1(3,3);
% Display the Fundamental matrix
disp('Fundamental matrix F:');
disp(F1);

%% STEP 5
V(:,1) = [100;-400;2000];
V(:,2) = [300;-400;3000];
V(:,3) = [500;-400;4000];
V(:,4) = [700;-400;2000];
V(:,5) = [900;-400;3000];
V(:,6) = [100;-50;4000];
V(:,7) = [300;-50;2000];
V(:,8) = [500;-50;3000];
V(:,9) = [700;-50;4000];
V(:,10) = [900;-50;2000];
V(:,11) = [100;50;3000];
V(:,12) = [300;50;4000];
V(:,13) = [500;50;2000];
V(:,14) = [700;50;3000];
V(:,15) = [900;50;4000];
V(:,16) = [100;400;2000];
V(:,17) = [300;400;3000];
V(:,18) = [500;400;4000];
V(:,19) = [700;400;2000];
V(:,20) = [900;400;3000];

%% STEP 6
% Projection on image planes
cam1_p2d = mvg_projectPointToImagePlane(V,P1);
cam2_p2d = mvg_projectPointToImagePlane(V,P2);

%% STEP 7
% example of the plotting functions 
% Draw 2D projections on image planes
cam1_fig = mvg_show_projected_points(cam1_p2d(1:2,:),imageSize,'Projected points on image plane 1');
cam2_fig = mvg_show_projected_points(cam2_p2d(1:2,:),imageSize,'Projected points on image plane 2');

%% step 8

[~,numcol] = size(V);

% Create U_new using a more compact matrix operation
U_new = [cam2_p2d(1,:)'.*cam1_p2d(1,:)' cam2_p2d(1,:)'.*cam1_p2d(2,:)' ...
    cam2_p2d(1,:)' cam2_p2d(2,:)'.*cam1_p2d(1,:)' cam2_p2d(2,:)'.*cam1_p2d(2,:)' ...
    cam2_p2d(2,:)' cam1_p2d(1,:)' cam1_p2d(2,:)' ones(numcol,1)];


% Compute SVD of U_new
[~, ~, V_eight] = svd(U_new);

% Reshape the last column of V_eight to form F1
F= reshape(V_eight(:,9), 3,3)';
% Compute SVD of F1
[U_one, D_one, V_one] = svd(F);
% Enforce rank-2 constraint on F1
F = U_one*diag([D_one(1,1) D_one(2,2) 0])*-V_one';

%% step 9
error = 0;
for i = 1:size(F1,1)
    for j = 1:size(F1,2)
        error = error + abs(F1(i,j) - F(i,j));
    end
end
% Display the difference
disp(['Sum of Absolute Differences: ', num2str(error)]);
%% step 10
% Draw epipolar lines
%[~,~,c1_l_coeff,c2_l_coeff] = mvg_compute_epipolar_geom_modif(cam1_p2d,cam2_p2d,F);
%[cam1_fig,cam2_fig] = mvg_show_epipolar_lines(cam1_fig, cam2_fig, c1_l_coeff,c2_l_coeff, [-400,1;300,400],'b');

% Draw epipoles
%[U,~,V] = svd(F);
%ep_1 = V(:,end); % The epipole in the first image is the last column of V
%ep_2 = U(:,end); % The epipole in the second image is the last column of U

% Normalize the epipoles (to make the third coordinate 1)
%ep_1 = ep_1 / ep_1(3);
%ep_2 = ep_2 / ep_2(3);

%[~,~] = mvg_show_epipoles(cam1_fig, cam2_fig,ep_1,ep_2);
%% step 11
proj_size = size(cam1_p2d);
r_one = normrnd(0,0.5,proj_size);
r_two = normrnd(0,0.5,proj_size);
cam1_p2d_noise = cam1_p2d + r_one;
cam2_p2d_noise = cam2_p2d + r_two;
%% step 12

% numPoints = size(V, 2);
% 
% % Create U_new using a more compact matrix operation
% U_new = [cam2_p2d_noise(1,:)'.*cam1_p2d_noise(1,:)', cam2_p2d_noise(1,:)'.*cam1_p2d_noise(2,:)', ...
%          cam2_p2d_noise(1,:)', cam2_p2d_noise(2,:)'.*cam1_p2d_noise(1,:)', cam2_p2d_noise(2,:)'.*cam1_p2d_noise(2,:)', ...
%          cam2_p2d_noise(2,:)', cam1_p2d_noise(1,:)', cam1_p2d_noise(2,:)', ones(numPoints, 1)];
% 
% % Compute SVD of U_new
% [~, ~, V_eight] = svd(U_new);
% 
% % Reshape the last column of V_eight to form F1
% F1_Noise = reshape(V_eight(:, end), 3, 3)';
% 
% % Compute SVD of F1
% [U_one, D_one, V_one] = svd(F1_Noise);
% 
% % Enforce rank-2 constraint on F1
% F1_Noise = U_one * diag([D_one(1,1), D_one(2,2), 0]) * V_one';

% error = 0;
% for i = 1:size(F1,1)
%     for j = 1:size(F1,2)
%         error = error + abs(F1(i,j) - F1_Noise(i,j));
%     end
% end
% % Display the difference
% disp(['Sum of Absolute Differences: ', num2str(error)]);
% 

% [~,~,c1_l_coeff,c2_l_coeff] = mvg_compute_epipolar_geom_modif(cam1_p2d,cam2_p2d,F1_Noise);
% [cam1_fig,cam2_fig] = mvg_show_epipolar_lines(cam1_fig, cam2_fig, c1_l_coeff,c2_l_coeff, [-400,1;300,400],'b');
% 
% [U,~,V] = svd(F1_Noise);
% ep_1 = V(:,end); % The epipole in the first image is the last column of V
% ep_2 = U(:,end); % The epipole in the second image is the last column of U
% 
% ep_1 = ep_1 / ep_1(3);
% ep_2 = ep_2 / ep_2(3);
% 
% [~,~] = mvg_show_epipoles(cam1_fig, cam2_fig,ep_1,ep_2);
%% step 13
proj_size = size(cam1_p2d);
r_one = normrnd(0,1,proj_size);
r_two = normrnd(0,1,proj_size);
cam1_p2d_noise = cam1_p2d + r_one;
cam2_p2d_noise = cam2_p2d + r_two;
[~,numcol] = size(V);
U_new = [cam2_p2d_noise(1,:)'.*cam1_p2d_noise(1,:)' cam2_p2d_noise(1,:)'.*cam1_p2d_noise(2,:)' cam2_p2d_noise(1,:)' cam2_p2d_noise(2,:)'.*cam1_p2d_noise(1,:)' cam2_p2d_noise(2,:)'.*cam1_p2d_noise(2,:)' cam2_p2d_noise(2,:)' cam1_p2d_noise(1,:)' cam1_p2d_noise(2,:)' ones(numcol,1)];

[~, ~,V_eight_noise] = svd(U_new);
F_eight_noise = reshape(V_eight_noise(:,9),3,3)';
[U_one_noise, D_one_noise, V_one_noise] = svd(F_eight_noise);
F_eight_noise=U_one_noise*diag([D_one_noise(1,1) D_one_noise(2,2) 0])*V_one_noise';
% 
% error = 0;
% for i = 1:size(F1,1)
%     for j = 1:size(F1,2)
%         error = error + abs(F1(i,j) - F_eight_noise(i,j));
%     end
% end
% 
% disp(['Sum of Absolute Differences: ', num2str(error)]);
% 
% [~,~,c1_l_coeff,c2_l_coeff] = mvg_compute_epipolar_geom_modif(cam1_p2d,cam2_p2d,F_eight_noise);
% [cam1_fig,cam2_fig] = mvg_show_epipolar_lines(cam1_fig, cam2_fig, c1_l_coeff,c2_l_coeff, [-400,1;300,400],'b');
% 
% 
% [U,~,V] = svd(F_eight_noise);
% ep_1 = V(:,end); % The epipole in the first image is the last column of V
% ep_2 = U(:,end); % The epipole in the second image is the last column of U
% 
% 
% ep_1 = ep_1 / ep_1(3);
% ep_2 = ep_2 / ep_2(3);
% 
% [~,~] = mvg_show_epipoles(cam1_fig, cam2_fig,ep_1,ep_2);
%% step 14


% Normalize points from camera 1
mean_cam1_p2d = mean(cam1_p2d(1:2,:), 2);
std_cam1_p2d = std(cam1_p2d(1:2,:), 0, 2);
T1 = [sqrt(2)/std_cam1_p2d(1) 0 -sqrt(2)/std_cam1_p2d(1)*mean_cam1_p2d(1); 
      0 sqrt(2)/std_cam1_p2d(2) -sqrt(2)/std_cam1_p2d(2)*mean_cam1_p2d(2); 
      0 0 1];
cam1_p2d_normalized = T1 * [cam1_p2d(1:2,:); ones(1, size(cam1_p2d, 2))];

% Normalize points from camera 2
mean_cam2_p2d = mean(cam2_p2d(1:2,:), 2);
std_cam2_p2d = std(cam2_p2d(1:2,:), 0, 2);
T2 = [sqrt(2)/std_cam2_p2d(1) 0 -sqrt(2)/std_cam2_p2d(1)*mean_cam2_p2d(1); 
      0 sqrt(2)/std_cam2_p2d(2) -sqrt(2)/std_cam2_p2d(2)*mean_cam2_p2d(2); 
      0 0 1];
cam2_p2d_normalized = T2 * [cam2_p2d(1:2,:); ones(1, size(cam2_p2d, 2))];

% Ensure that points are in homogeneous coordinates
cam1_p2d_normalized = cam1_p2d_normalized ./ cam1_p2d_normalized(3,:);
cam2_p2d_normalized = cam2_p2d_normalized ./ cam2_p2d_normalized(3,:);

%for step 8:
U_new_normalized = [cam2_p2d_normalized(1,:)'.*cam1_p2d_normalized(1,:)' cam2_p2d_normalized(1,:)'.*cam1_p2d_normalized(2,:)' ...
    cam2_p2d_normalized(1,:)' cam2_p2d_normalized(2,:)'.*cam1_p2d_normalized(1,:)' cam2_p2d_normalized(2,:)'.*cam1_p2d_normalized(2,:)' ...
    cam2_p2d_normalized(2,:)' cam1_p2d_normalized(1,:)' cam1_p2d_normalized(2,:)' ones(numcol,1)];

% Compute SVD of U_new_normalized
[~, ~, V_eight_normalized] = svd(U_new_normalized);

% Reshape the last column of V_eight_normalized to form F_normalized
F_normalized = reshape(V_eight_normalized(:,9), 3,3)';

% Enforce the rank-2 constraint on F_normalized
[U_normalized, D_normalized, V_normalized] = svd(F_normalized);
F_normalized = U_normalized * diag([D_normalized(1,1) D_normalized(2,2) 0]) * V_normalized';

% Denormalize the F matrix
F_D = T2' * F_normalized * T1;

% Check the condition number of U_new_normalized
cond_normalized = cond(U_new_normalized);

% For comparison, compute the condition number of the original U_new
cond_unnormalized = cond(U_new);

% Display the condition numbers
disp(['Condition number with normalization: ', num2str(cond_normalized)]);
disp(['Condition number without normalization: ', num2str(cond_unnormalized)]);

error = 0;
for i = 1:size(F_eight_noise,1)
    for j = 1:size(F_eight_noise,2)
        error = error + abs(F_eight_noise(i,j) - F_D(i,j));
    end
end
disp(['Sum of Absolute Differences: ', num2str(error)]);

[~,~,c1_l_coeff,c2_l_coeff] = mvg_compute_epipolar_geom_modif(cam1_p2d,cam2_p2d,F_D);
[cam1_fig,cam2_fig] = mvg_show_epipolar_lines(cam1_fig, cam2_fig, c1_l_coeff,c2_l_coeff, [-400,1;300,400],'b');

[U,~,V] = svd(F_D);
ep_1 = V(:,end); % The epipole in the first image is the last column of V
ep_2 = U(:,end); % The epipole in the second image is the last column of U

ep_1 = ep_1 / ep_1(3);
ep_2 = ep_2 / ep_2(3);

[~,~] = mvg_show_epipoles(cam1_fig, cam2_fig,ep_1,ep_2);

return;
