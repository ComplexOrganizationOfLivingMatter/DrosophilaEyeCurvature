function coordinatesEyeSurface = ObtainAdultEyeFromVolume(coordinates,boundaryPoints,S)

%% Inputs

% coordinates: 3D coordinate matrix that represents the volume of 
%              adultEyeSegmentationImage in microns interpolated.

% boundaryPoints: Vector with the index of the points of coordinates that
%                 belong to its boundary in the xy plane.

% S: alphaShape object that contains the information of the triangulation
%    of the volume matrix of coordinates.

%% Output

% coordinatesEyeSurface: 3D coordinate matrix that represents the 3D
%                        surface of the adult eye.

%% Code explanation

% This code will obtain the 3D surface of the adult eye 
% (coordinatesEyeSurface) from its volume coordinate matrix (coordinates).

% In order to do that, we will identify the apex, which is the point of
% highest z-coordinate. Then, we will identify the points connected to that
% points by the triangles of the mesh. All these points will be identified
% as valid points. We will apply the algorithm iterativetly until we reach
% the boundaryPointsm that will work as a stopping line, since they will
% stop the algorithm (if a boundaryPoint is reached, we don't keep on
% identifying valid points. Finally, the valid points and boundaryPoints
% will be added to coordinatesEyeSurface.

%% Identifying the apex and the boundary points in the mesh

% First, we obtain the highest z-coordinate point and the connectivity 
% matrix (elements) of the mesh, that will be the boundaryFacets of S:

rowIndexApex = find(coordinates(:,3) == max(coordinates(:,3)));
elements = S.boundaryFacets;

% Now, we will identify the points in the mesh boundaryPoints belong to:

rowIndexBoundaryPoints = [];

for i = 1:length(boundaryPoints)

    [rowIndexBoundaryPointsIteration,~] = find(elements == boundaryPoints(i));
    rowIndexBoundaryPoints = [rowIndexBoundaryPoints;rowIndexBoundaryPointsIteration];

end

boundaryPointsMesh = unique(elements(rowIndexBoundaryPoints,:));

%% Applyin the algorithm to extract the adult eye surface

% To implement the algorithm we need to define two matrices. The first is
% externalSurface (since we will obtain, from the pov of the apex, the
% external surface, not the interior) and the externalSurfaceAuxiliar. A
% detailed information of both follows in the next lines:

% externalSurface: matrix that contains the points of the eye.

% externalSurfaceAuxiliar : matrix that contains the row of points whose
%                           neighbours we want to study. We will empty this
%                           matrix in every interation.

% So, keeping in mind that definitions, we initialize them as:

externalSurface = row_z_highest;
externalSurfaceIndexAuxiliar = row_z_highest;

% Also, we initialize the points that are connected by the triangles in the
% mesh with externalSurfaceAuxiliar:

rowsIndexSharePoints = [];

% It's important to remark that the algorithm will apply iteratively until
% the condition (there is no more points to identify) is reached, so we set
% an iterationsNumber high enough without any additional computational 
% cost:

iterationsNumber = 1e4;

for i = 1:iterationsNumber

    % First, we need to find the triangles that share the points in 
    % externalSurfaceIndexAuxiliar

    for j = 1:size(externalSurfaceIndexAuxiliar,1)

        [rowsIndexSharePointsIteration,~] = find(elements == externalSurfaceIndexAuxiliar(j));
        rowsIndexSharePoints = [rowsIndexSharePoints;rowsIndexSharePointsIteration];

    end

    % Using unique we obtain the points that belong to the triangles shared
    % by the points in externalSurfaceIndexAuxiliar:

    studyPoints = unique(elements(rowsIndexSharePoints,:)); 

    % From this, we delete all the points that we already know they satisfy
    % our conditons (in the first interation for example it will be the
    % point with higher z):

    pointsAlreadyStudied = [];

    for j = 1:size(studyPoints,1)

        if ~isempty(find(externalSurface == studyPoints(j)))

            pointsAlreadyStudied = [pointsAlreadyStudied;j];

        end

    end

    studyPoints(pointsAlreadyStudied) = [];    

    % Now, we check if some of the study points are points that belong to
    % the boundary and delete them:

    pointsBelongBoundary = [];

    for j = 1:size(studyPoints,1)

        if ~isempty(find(boundaryPointsMesh == studyPoints(j)))

            pointsBelongBoundary = [pointsBelongBoundary;j];

        end

    end
    
    studyPoints(pointsBelongBoundary) = [];

    % Now, we impose the break condition, which is that the matrix of
    % studyPoints in empty (we have deleted from it the points that have
    % been studied previously and the points in the boundary, so in some
    % moment, there won't be more points to study):

    if isempty(studyPoints)

        break

    else

        externalSurface = [externalSurface;studyPoints];
        externalSurfaceIndexAuxiliar = studyPoints;

    end 

end

%% Obtaining the externalSurface and coordinatesEyeSurface

% Finally, we extract the coordinates from externalSurface and add the 
% ones in the boundary:

externalSurface = [externalSurface;boundaryPoints];
coordinatesEyeSurface = coordinates(ext_surf,:);