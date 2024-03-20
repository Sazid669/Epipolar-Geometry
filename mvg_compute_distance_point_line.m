function [ res_d, d  ] = mvg_compute_distance_point_line( lines,points )
%COMPUTE_DISTANCE_POINT_LINE Compute statistic of distance between point
%and corresponding lines
d = zeros(size(lines,2),1); %preallocating
for ii=1:size(lines,2)
   d(ii) = abs(lines(:,ii)'*points(:,ii))/sqrt(lines(1,ii)*lines(1,ii)+lines(2,ii)*lines(2,ii));  
end
    res_d(1) = sum(d);
    res_d(2) = mean(d);
    res_d(3) = std(d);
end

