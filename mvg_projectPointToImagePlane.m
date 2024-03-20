function [ p_2d ] = mvg_projectPointToImagePlane( p_3d, P )
%PROJECTPOINT project point onto image plane using projection matrix P
if size(p_3d,1) == 3
    p_3d = [p_3d; ones(1,size(p_3d,2))];
end
    
    p_2d_h = P * p_3d;
    p_2d(1,:) = p_2d_h(1,:)./p_2d_h(3,:);
    p_2d(2,:) = p_2d_h(2,:)./p_2d_h(3,:);

end

