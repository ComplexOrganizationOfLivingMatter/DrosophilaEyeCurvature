%% Code explanation

% This code will build the uniform identical triangles mesh version of sWT
% simulation. As we did for sWT simulation, this mesh will accomplish the
% requirements of the computational model, e.g, being a quadratic
% triangulation, define the Dirichlet Boundary Conditions and the Nodal
% points.

%%  Material Properties

% Wed define the porperties of the isometric material we will simulate and
% the pressure apllied to it;

young = 6.825e9;
poiss = 3e-1;
denss = 0;
thick = 0.125;
p = 1e5;

%% Loading necessary data

% To create an uniform identical triangles pattern version of the sWT
% pattern we need to load the data from the sWT simulation to set the
% frontier inside wich we will generate an uniform lattice:

load("sWT_1_Male_BuiltMesh.mat",'k','coordinates','conversionFactor','elements');

% To have a good first guess of the number of triangles that the uniform
% mesh will have (the same as sWT ideally). We need to load tge data from
% the segmentation in order to have an area idea:

load("sWT_1_Male_MeshVariables.mat",'validTriangles','inverseTriangles')

%% Computing the initial dimensions of uniform triangles

% We could apply the algorithm (see below) to get the size of the triangles
% with an initial random guess in the base length and height of triangles.
% However, it will increase so much the computational cost.

% For this reason, we will make a first guess of the dimensions of the
% identical triangles, calculating the total area of the region inside the
% boundary of sWT and then dividing by the total number of triangles
% (validTriangles + inverseTriangles) to determine the initial guess of
% area.

% So, the area of WT, the total number of triangles and the initial guess
% of area are:

sWTArea = polyarea(coordinates(k,1),coordinates(k,2));
trianglesNumber = size(validTriangles,1) + size(inverseTriangles,1);
trianglesArea = sWTArea/trianglesNumber;

% Now, we obtain the height and base of the identical triangles:

[baseIdenticalTriangles,heightIdenticalTriangles] = ObtainBaseHeightIdenticalTriangles(validTrianglesCoordinates,validTriangles,triangleArea,conversionFactor);

%% Generating the uniform lattice and the uniform mesh

% Once we have the inital guesses for the base and height we could obtain
% the uniform lattice and the uniform mesh:

[uniformTrianglesCoordinates,uniformTriangles,P] = GetUniformMesh(baseIdenticalTriangles,heightIdenticalTriangles,coordinates,elements,k);

%% Building the mesh for the computational model

% We will follow the same process that we did for the sWT.

[coordinates,elements] = FromMeshSegmentationToQuadraticTriangulation(uniformTriangles,uniformTrianglesCoordinates,P,conversionFactor);

%% Dirichlet Boundary Conditions, curations and nodal points

% To obtain the Dirichlet Boundary Conditions we need to identify the
% boundary points of the mesh, excluding the midpoints of the triangles
% that belong to the boundary:

coordinatesVertices = coordinates(1:size(P,1),1:2);
k = boundary(coord_vert,1); % be careful with the number

% In this case, due to the process of uniform lattice, it's possible to
% have a point in coordinates that don't appear in elements (a point which
% is alone in the boundary and impossible to triangulate). So, we need to 
% delete the point and fix the elements matrix:

pointsToDelete = setdiff(k,elements);
pointsToDelete = sort(pointsToDelete);
coordinates(points_delete,:) = [];
elements(elements>points_delete(1)) = elements(elements>points_delete(1)) - 1;

% Once it has been corrected, we need to recalculate the Dirichlet Boundary
% Conditions:

elementsVertices = unique(elements(:,[1,3,5]));
coordinatesVertices = coordinates(elementsVertices,1:2);
k = boundary(coordinatesVertices,0.99);

% Once the boundary has been identified, the Dirichlet Boundary Conditions
% are defined in fixnodes and the nodal points in pointload:

fixnodes = [k,ones(size(k,1),1),zeros(size(k,1),1);k,2*ones(size(k,1),1),zeros(size(k,1),1);k,3*ones(size(k,1),1),zeros(size(k,1),1)];
pointloadnodes = setdiff(elementsVertices,k);
pointload = [pointloadnodes,3*ones(size(pointloadnodes,1),1),p*ones(size(pointloadnodes,1),1)];

%% Orient the mesh

% As sWT was previously oriented, we don't need orientation for sIT.

% Side loads

uniload = sparse(size(elements,1),1);

% Saving the mesh

save("sIT_1_Male_BuiltMesh.mat")