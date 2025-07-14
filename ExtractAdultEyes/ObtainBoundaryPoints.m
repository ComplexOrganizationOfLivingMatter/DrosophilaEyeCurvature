function boundaryPoints = ObtainBoundaryPoints(coordinates)

%% Inputs

% coordinates: 3D coordinate matrix that represents the volume of 
%              adultEyeSegmentationImage in microns interpolated.

%% Output

% boundaryPoints: Vector with the index of the points of coordinates that
%                 belong to its boundary in the xy plane.

%% Code explanation

% This code will obtain the index of the boundary points of coordinates
% matrix in the xy plane.

%% Code

% To obtain the boundary we use the boundary function:

boundaryPoints = boundary(coordinates(:,1),coordinates(:,2),1);

% In this case the first and the last element of boundaryPoints are the
% because its a closed curve so we keep the unique values:

boundaryPoints = unique(boundaryPoints);

