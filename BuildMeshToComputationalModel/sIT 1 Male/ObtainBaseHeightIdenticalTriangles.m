function [baseIdenticalTriangles,heightIdenticalTriangles] = ObtainBaseHeightIdenticalTriangles(validTrianglesCoordinates,validTriangles,triangleArea,conversionFactor)

%% Inputs

% validTrianglesCoordinates: Coordinate tensor wich contains the
%                            coordinates of the validTriangles. It's a
%                            tensor because it makes easier the
%                            postprocessing.

% validTriangles: Connectivity matrix that represents the triangles that
%                 correspond to the set of valid triangles.

% triangleArea: Area of the uniform triangles considering just total area
%               of sWT and the total number of triangles in the
%               segmentation data.

% conversionFactor: It's the conversion factor from pixels² to microns².
%                   It's the square power of the spatial revolution of 
%                   basalImage.

%% Outputs

% baseIdenticalTriangles: The initial guess of the base length of the
%                         uniform triangles.

% heightIdenticalTriangles: The initial guess of the height of the uniform
%                           triangles.

%% Code explanation

% This code will obtain the estimation of the base and height of identical
% triangles based of validTrianglesCoordinates and its connectivity matrix
% validTriangles.

%% Transform the coordinates to microns and obtain the base and height:

validTrianglesCoordinates = validTrianglesCoordinates*sqrt(conversionFactor);

for i = 1:size(validTriangles,1)

    bases(i) = sqrt((validTrianglesCoordinates(i,1,1) - validTrianglesCoordinates(i,2,1))^2 + (validTrianglesCoordinates(i,1,2)- validTrianglesCoordinates(i,2,2))^2);
    aux_1 = abs((validTrianglesCoordinates(i,2,1) - validTrianglesCoordinates(i,1,1))*(validTrianglesCoordinates(i,1,2) - validTrianglesCoordinates(i,3,2)) - (validTrianglesCoordinates(i,1,1) - validTrianglesCoordinates(i,3,1))*(validTrianglesCoordinates(i,2,2) - validTrianglesCoordinates(i,1,2)));
    aux_2 = sqrt((validTrianglesCoordinates(i,2,1) - validTrianglesCoordinates(i,1,1))^2 + (validTrianglesCoordinates(i,2,2) - validTrianglesCoordinates(i,1,2))^2);
    heights(i) = aux_1/aux_2;

end

% So, the mean base, mean height and mean area are:

meanBase = mean(bases);
meanHeight = mean(heights);
meanArea = meanBase*meanHeight/2;

%% Correction to meanBase

% Empiracally, its shown that exists a difference between triangleArea and
% meanArea it's necessary to correct meanBase for not losing triangles
% since if meanArea is bigger than triangleArea some triangles will be
% lost.

correction = 2*(triangleArea-meanArea)/meanHeight;
baseIdenticalTriangles = meanBase + correction;
heightIdenticalTriangles = meanHeight;

