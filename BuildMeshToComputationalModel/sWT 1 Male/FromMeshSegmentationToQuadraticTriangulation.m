function [coordinates,elements] = FromMeshSegmentationToQuadraticTriangulation(validTriangles,validTrianglesCoordinates,P,conversionFactor)

%% Inputs

% validTriangles: Connectivity matrix that represents the triangles that
%                 correspond to the set of valid triangles.

% validTrianglesCoordinates: Coordinate tensor wich contains the
%                            coordinates of the validTriangles. It's a
%                            tensor because it makes easier the
%                            postprocessing.

% P: Coordinate matrix corresponding to the set of valid triangles and the
%    inverse triangles.

% conversionFactor: It's the conversion factor from pixels² to microns².
%                   It's the square power of the spatial revolution of 
%                   basalImage.

%% Outputs

% coordinates: Coordinate matrix with the coordinates of the points that
%              represent the quadratic triangulation of 
%              validTrianglesCoordinates

% elements: Connectivity matrix of coordinates.

%% Code explanation

% This code will obtain the coordinates matrix and its connectivity matrix,
% elements, to transform the validTrianglesCoordinates in a quadratic
% triangulation.

%% Obtaining the middle points and coordinates

% We concatenate the first column of validTrianglesCoordinates to
% validTrianglesCoordinates. In that way, the third colum will obtain the
% midpoint with the first on and the bucle will be easier. So, we will
% create validTrianglesCoordinates_aux:

validTrianglesCoordinates_aux = [validTrianglesCoordinates,validTrianglesCoordinates(:,1,:)];

for i = 1:size(validTrianglesCoordinates,1)

    for k = 1:3

        CoordinatesMidPoints(i,k,:) = (validTrianglesCoordinates_aux(i,k,:) + validTrianglesCoordinates_aux(i,k+1,:))/2;

    end

end

% Now we need the coordinates in a matrix, not in a tensor. We will have P
% and below we will have the coordinates of the midpoints where the first 3
% rows will be the 3 midpoints of the first triangle. So, we take the x and
% y coordinates and then we traspose:

CoordinatesMidPoints_X = CoordinatesMidPoints(:,:,1).';

CoordinatesMidPoints_Y = CoordinatesMidPoints(:,:,2).';

% Now we have a matrix where column i is the row i of CoordinatesMidPoints.
% Now we concatenate in a column vector all the colums.

coordinates_X = CoordinatesMidPoints_X(:,1,:);
coordinates_Y = CoordinatesMidPoints_Y(:,1,:);

% We transform the all the tensor in a matrix:

for i = 2:size(CoordinatesMidPoints_X,2)

    coordinates_X = [coordinates_X;CoordinatesMidPoints_X(:,i)];
    coordinates_Y = [coordinates_Y;CoordinatesMidPoints_Y(:,i)];

end

% Now we make the coordinate matrix of midpoints in 2D:

coordinatesMidPoints2D = [coordinates_X,coordinates_Y];

% To obtain the coordinates matrix, we add coordinatesMidPoints2D to P and
% transform the matrix to microns:

coordinates = [P;coordinatesMidPoints2D]*sqrt(conversionFactor);

% We finally transform the coordinates to 3D:

coordinates = [coordinates,zeros(size(coordinates,1),1)];

%% Obtaining the connectivity matrix elements

% Now we need the elements, that are the index of the vertex, in other
% words, the row in the coordinate matrix that represents the vertex or
% node of the triangle. In this way:

for i= 1:size(validTrianglesCoordinates,1)

    Index_Matrix(i,:) = [validTriangles(i,:),[3*(i-1)+1+size(P,1):3*(i-1)+3+size(P,1)]];

end

% In this way, we have the first 3 columns (1,2,3) that belongs to the vertices of
% the triangle and the 3 last (4,5,6) that belong to the midpoints. As we
% have calculated the midpoints in the way, the first midpoint is the one
% between vertex 1 and 2, in elements, the second row of indices should the
% one row for mid points, so:

elements = [Index_Matrix(:,1),Index_Matrix(:,4),Index_Matrix(:,2),Index_Matrix(:,5),Index_Matrix(:,3),Index_Matrix(:,6)];
