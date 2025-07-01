function [GaussianCurvaturePercentiles,Coordinates,points_in_percentil,K] = GetGaussianCurvaturePercentiles(Coordinates,K,n)

%% Inputs

% Coordinates: Matrix of points m x 3, with n being the number of points:

% K: Gaussian curvature of the points in Coordinates. It's vector of length
%    m.

% n: The number of regions in which de distribution of distance is divided.

%% Outputs

% GaussianCurvaturePercentiles: It's a vector of length n that contains the
%                               summation of K of points that belong to the
%                               distance percentiles.

%% Code

% We need to apply a threshold to K because there are points that haven't
% been well calculated. So we will obtain the mean of K and the std. All
% the points whose abs is greater than 0.2 will be deleted. Then that are above percentil 99.1 will be deleted:

% pointsToDeleteAbs = find(abs(K)>0.1);
% K(pointsToDeleteAbs) = [];
% Coordinates(pointsToDeleteAbs,:) = [];
% percentilThresholdAbove = prctile(K,99);
% pointsDeleteAbove = find(K >= percentilThresholdAbove);
% percentilThresholBelow = prctile(K,1);
% pointsDeleteBelow = find(K <= percentilThresholBelow);
% pointsDelete = [pointsDeleteAbove;pointsDeleteBelow];
% % pointsDelete = find(isoutlier(K));
% Coordinates(pointsDelete,:) = [];
% K(pointsDelete) = [];

% We will project the points into 2D plane in polar coordinates and will
% divide them in regions, percentiles, by the distance to the center of
% their distribution. Every percentil must have the same number of points.
% Finally, we will sum up the total Gaussian curvature of all points in
% every region to obtain T_G_C_b_P_P.

x = Coordinates(:,1);

y = Coordinates(:,2);

x_center = mean(Coordinates(:,1));

y_center = mean(Coordinates(:,2));

% Now we obtain distance r and the angle theta of each point:

for i = 1:length(x)
    r(i) = sqrt((x(i) - x_center)^2 + (y(i) - y_center)^2);
    distance_x_to_origin = x(i) - x_center;
    distance_y_to_origin = y(i) - y_center;
    theta(i) = atan(distance_y_to_origin/distance_x_to_origin);
end

% We obtain the percentiles of the distribution and the value of r that
% defines them. We will divide the percentiles by 100 to normalize them
% between [0,1].

j = 1:n;

r_percentiles = ((100*(j-1))/n + 1);

% Now, we obtain the values:

r_percentiles_values = prctile(r,r_percentiles);

% We need to identify the points that belongs to each percentil and then we
% will sum up the total gaussian curvature of this region to save it in a
% vector, in T_G_C_b_P_P:

for i = 1:n+1

    if i == 1

        points_in_percentil{i} = find(r<r_percentiles_values(i));
        GaussianCurvaturePercentiles(i) = sum(K(points_in_percentil{i}));

    else

        if i ~= n+1
            
            points_in_percentil{i} = find(r>=r_percentiles_values(i-1) & r<r_percentiles_values(i));
            GaussianCurvaturePercentiles(i) = sum(K(points_in_percentil{i}));

        else

            points_in_percentil{i} = find(r>=r_percentiles_values(i-1));
            GaussianCurvaturePercentiles(i) = sum(K(points_in_percentil{i}));

        end

    end

end
