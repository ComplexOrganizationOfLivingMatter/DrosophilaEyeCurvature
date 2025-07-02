function [uniformValidTrianglesCoordinates,uniformValidTriangles] = GetUniformTriangulation(uniformCoordinatesTensor,P,Tri)

%% Inputs

% uniformCoordinatesTensor: Coordinate tensor of Mx3x2. Each row 
%                           represents the 3 vertices of one triangle
%                           The x-y coordinates of each vertex are
%                           inside the 2 tensor dimension.

% P: Coordinate matrix of the uniform triangulation.

% Tri: Connectivity matrix that represent P.

%% Outputs

% uniformValidTrianglesCoordinates: Coordinate tensor of the
%                                   uniformaLattice set of validTriangles.

% uniformValidTriangles: Connectivity matrix that represents
%                        uniformValidTrianglesCoordinates.

%% Code explanation

% This code, based on the geometric features of the uniform lattice
% triangulation will obtain its set of validTriangles.

%% Algorithm to extract the uniformValidTriangles set

uniformValidTriangles = [];

for j = 1:size(Tri,1)

    maxYCoordinates = max(uniformCoordinatesTensor(j,:,2));
    comparedCoordinates = uniformCoordinatesTensor(j,:,2) == maxYCoordinates;

    if sum(comparedCoordinates) == 1

        uniformValidTriangles = [uniformValidTriangles;Tri(j,:)];

    end
end

%% Algorithm to extract the uniformValidTrianglesCoordinates

uniformValidTrianglesCoordinates = getCoordinatesTensor(P,uniformValidTriangles);