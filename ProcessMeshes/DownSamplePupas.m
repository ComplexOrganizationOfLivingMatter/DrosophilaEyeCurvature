function [ConnectivityMatrixDownsampled,CoordinatesDownsample] = DownSamplePupas(coordinates)

% This code will downsample the point cloud and the will triangulate the
% points:

ptCloud = pointCloud(coordinates);

% Calculate the number of points to sample using percentageFactor (a
% percentageFactor of 0.05 means that we keep just the 5% of the points.
% We will also determine the minPointsNumber which is the minimum
% number of points just in case we use a percentageFactor too low.

percentageFactor = 0.05;
minPointsNumber = 5000;
numOriginalPoints = size(ptCloud.Location, 1);
numPoints = max(minPointsNumber, round(numOriginalPoints * percentageFactor)); % Ensure at least 100 points

% Use farthest point sampling to get a uniform distribution of points

idx = farthestPointSampling(ptCloud, numPoints);

% We know obtain the points and triangulate them:

CoordinatesDownsample = ptCloud.Location(idx, :);
tri2D = delaunay(CoordinatesDownsample(:,1),CoordinatesDownsample(:,2));
dt = triangulation(tri2D,CoordinatesDownsample(:,1),CoordinatesDownsample(:,2),CoordinatesDownsample(:,3));
CoordinatesDownsample = dt.Points;
ConnectivityMatrixDownsampled = dt.ConnectivityList;
