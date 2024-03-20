function [ h ] = mvg_show_projected_points( points_2d,image_size,titleStr )
% [ h ] = mvg_show_projected_points( points_2d,image_size,titleStr )
%
% Show projected points on image plane together with image borders. 
% Returns a handle to the figure created, so that other functions can draw
% on top

%some may be outside of image - find factor that encapsulates all points
factors(1,:) = points_2d(1,:) ./ image_size(1);
factors(2,:) = points_2d(2,:) ./ image_size(2);

%no matter what we show the whole image [0,1] factors
max_f = max(1,max(factors(:)))-1; % so the two factors are comparable
min_f = min(0,min(factors(:)));
%find the factor that encapsulates all points and scale a bit more
factor = 1.25*max(abs(min_f),max_f);

h = figure;
axis([-factor*image_size(1) (1+factor)*image_size(1) -factor*image_size(2) (1+factor)*image_size(2)]);
axis ij;

title(titleStr);
hold on;
%show the border of image
rectangle('Position',[0,0,image_size(1),image_size(2)]);   
%show projected points
scatter(points_2d(1,:),points_2d(2,:),30,'m','x','LineWidth',1.5);
scatter(points_2d(1,:),points_2d(2,:),30,'m','o');
hold off;

end

