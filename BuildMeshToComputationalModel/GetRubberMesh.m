function [coordinates,elements] = GetRubberMesh(coordinatessWT,k)

%% Inputs

% coordinatessWT: Coordinate matrix that represents the sWT mesh.

% k: Vector whose components are the indexes of the points of 
%    coordinatessWT that belong to it's boundary.

%% Outputs

% coordinates: Coordinate matrix that represents the rubber-like mesh.

% elements: Connectivity matrix that represents the coordinates matrix.

%% Code explanation

% This code will obtain the coordinates and elements matrices satisfying
% the computational model requirements.

%% Algorithm to obtain rubber-like model

% We first obtain the x-y plane coordiantes of coordinatessWT and the we
% build the model:

coordinatessWT = coordinatessWT(:,1:2);
model = createpde;

% Now, we need to create the geometry with the boundary points. To do that,
% we need to create the geometry for a polygon, and this implies that the
% first row needs to be a 2. The next row os the number of line segments in
% the boundary of the polygon, which is length(coordinates(k,1)). The next
% n rows (n being the number of points) are the x-coordinates of the points
% in the boundary, and the next n rows are the y-coordinates of the points
% in the bounary.

boundary_matrix = [2;length(coordinatessWT(k,1));coordinatessWT(k,1);coordinatessWT(k,2)];
gm = boundary_matrix;

% Now, we need to build the descomposed geometry matrix g:

sf = 'boundary_matrix';
ns = char('boundary_matrix')';
g = decsg(gm,sf,ns);

% Using this matrix we can create the geometry and the mesh:

geometry = geometryFromEdges(model,g);
m1 = generateMesh(model,"Hmax",10);

%% Extract the coordinates and elements matrices

% From m1 we can easily obtain these matrices:

coordinates = m1.Nodes;
coordinates = [coordinates,zeros(size(coordinates,1),1)];
elements = m1.Elements';
elements = elements(:,[1,4,2,5,3,6]);

