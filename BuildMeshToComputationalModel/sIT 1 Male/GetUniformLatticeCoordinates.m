function uniformLatticeCoordinates = GetUniformLatticeCoordinates(basePoints,basePointsDisplaced,heightPoints)

%% Inputs

% basePoints: Vector that have the x-coordinates of the uniform lattice odd
%             rows.

% basePointsDisplaced: Vector that have the x-coordinates of the uniform
%                      lattice even rows.

% heightPoints: Vector that have the y-coordinates of all rows of the
%               uniform lattice.

%% Output

% uniformLatticeCoordinates: Coordinate matrix that represents the uniform 
%                            lattice.

%% Code explanation

% This code will obtain the coordinate matrix of the uniform lattice using
% the distribution of points of basePoints, basePointsDisplaced and 
% heightPoints.

uniformLatticeCoordinates = [];

for i = 1:length(heightPoints)/2

    firstRow = [basePoints',heightPoints(2*i-1)*ones(size(basePoints',1),1)];
    secondRow = [basePointsDisplaced',heightPoints(2*i)*ones(size(basePointsDisplaced',1),1)];
    uniformLatticeCoordinates=[uniformLatticeCoordinates;firstRow;secondRow];
    
end



