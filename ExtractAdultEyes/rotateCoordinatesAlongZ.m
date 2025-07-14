function rotatedCoordinates = rotateCoordinatesAlongZ(coordinates, direction)

%% Input
 
% coordinates: 3D coordinate matrix that represent the volume of 
%              adultEyeSegmentationImage in microns.

% direction: string, being 'x' or 'y'. The second highest eigenvector of
%            coordinates will be align along 'x' or 'y' direction.

%% Output

% rotatedCoordinates: 3D coordinate matrix oriented.

%% Code explanation

% This code rotates the point cloud of the highest eigenvector to be along
% z-direction. For the second highest eigenvector, it's possible to choose
% 'x' or 'y' direction.

%% Bring the point cloud center to the origin

coordinates = coordinates - mean(coordinates);

%% Align u normal vector along z-direction

% First, we obtain the highest eigenvalue and calculate the angles of the
% normal vector:

u = pcaEig(coordinates, 'max');
[alpha, beta] = unitVectorToAngle(u);       

% Now, we align the point cloud along x-axis followed by aligning along
% z-axis. To do that we need the rotational matrix:

[~, Ry, Rz] = rotationalMatrix(-alpha, pi-beta);
rotatedCoordinates = rotatePC(coordinates, Ry, Rz);

%% Align v (2nd highest eigenvector) normal vector along x or y direction

% We now decide the offset in function of direction:

switch direction
    case 'x'
        offset = 0;
    case 'y'
        offset = pi/2;
end 

% Now, we obtain the eigenvector of the second highest eigenvalue and
% calculate the projected v-vector along the xy-plane with respect to the
% x-axis:

v = pcaEig(rotatedCoordinates, 'middle');
[alpha, ~] = unitVectorToAngle(v);
[~, Ry, Rz] = rotationalMatrix(offset - alpha, 0);

%% Final rotation of the pointcloud

rotatedCoordinates = rotatePC(rotatedCoordinates, Ry, Rz);