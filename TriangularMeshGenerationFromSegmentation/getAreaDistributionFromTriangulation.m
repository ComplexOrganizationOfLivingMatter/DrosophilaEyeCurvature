function areaDistribution = getAreaDistributionFromTriangulation(trianglesCoordinates,conversionFactor)

%% Input

% trianglesCoordinates: Coordinate tensor wich contains the coordinates of
%                       the set of triangles that it represents. It's a
%                       Mx3x2 matrix, being M the number of triangles, 3
%                       because a triangle has 3 vertices and 2 because
%                       each vertex has a x-coordinate and an y-coordinate.

% conversionFactor: Conversion factor between the coordinates in pixels and
%                   the coordinates in microns of trianglesCoordinates.

%% Output

% areaDistribution: Area distribution of trianglesCoordinates. 

%% Code

areaDistribution = zeros(size(trianglesCoordinates,1),1);

for i=1:size(trianglesCoordinates,1)

    pgon = polyshape([trianglesCoordinates(i,:,1)],[trianglesCoordinates(i,:,2)]);
    areaDistribution(i) = area(pgon)*conversionFactor;

end