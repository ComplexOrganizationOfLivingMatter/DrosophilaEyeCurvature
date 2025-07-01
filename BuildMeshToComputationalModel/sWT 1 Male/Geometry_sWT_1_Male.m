%% Code Explanation

% This code will build the triangular mesh in order to accomplish the
% requierements of the computational model, e.g, being a quadratic
% triangulation, define the Dirichlet Boundary Conditions and the Nodal
% points.

% As every triangulation need manual curations since there is not a general
% way of identifying the Dirichlet Boundary Conditions and also the data
% from the segmentation involves artifacts. As a consequence, every mesh
% will have its own code and it won't exist a function to build the meshes
% for the computational model.

%%  Material Properties

% Wed define the porperties of the isometric material we will simulate and
% the pressure apllied to it;

young = 6.825e9;
poiss = 3e-1;
denss = 0;
thick = 0.125;
p = 1e5;

%% Creating the quadratic triangulation

% The computational model needs a quadratic triangulation (nodes in the 
% midlle point of edges of the triangles). To do that, we load the
% triangular mesh from the segmentations.

load('sWT_1_Male_MeshVariables.mat','validTriangles','validTrianglesCoordinates','P','conversionFactor')

% Now, from the tensor of validTrianglesCoordinates, we will obtain the
% coordinates matrix and the elements, both in a quadratic triangulation
% form:

[coordinates,elements] = FromMeshSegmentationToQuadraticTriangulation(validTriangles,validTrianglesCoordinates,P,conversionFactor);

%% Dirichlet Boundary Conditions and nodal points

% To obtain the Dirichlet Boundary Conditions we need to identify the
% boundary points of the mesh, excluding the midpoints of the triangles
% that belong to the boundary:

coordinatesVertices = coordinates(1:size(P,1),1:2);
k = boundary(coord_vert,0.75); % be careful with the number

% Once the boundary has been identified, the Dirichlet Boundary Conditions
% are defined in fixnodes and the nodal points in pointload:

fixnodes = [k,ones(size(k,1),1),zeros(size(k,1),1);k,2*ones(size(k,1),1),zeros(size(k,1),1);k,3*ones(size(k,1),1),zeros(size(k,1),1)];
pointloadnodes = setdiff(unique(validTriangles),k);
pointload = [pointloadnodes,3*ones(size(pointloadnodes,1),1),p*ones(size(pointloadnodes,1),1)];

%% Orient the mesh (Optional and manual since blinding in image adquisition)


% In order to have the mesh in the common orientation (D-V.A-P) we
% transform the coordinates matrix to orient it in the same way:

eyeCloud = pointCloud(coordinates);
thetha = 90; %counter-clockwise
rot = rotz(thetha);
trans = [0,0,0];
tform = rigid3d(rot,trans);
eyeCloudRotated = pctransform(eyeCloud,tform);
coordinates = eyeCloudRotated.Location;


% Now, we invert the coordinates:

coordinates2 = coordinates;
yCentroid = mean(coordinates(:,2));
coordinates2(:,2) = coordinates2(:,2) - yCentroid;
coordinates2(:,2) = -coordinates2(:,2);
coordinates = coordinates2;

% Side loads

uniload = sparse(size(elements,1),1);

%% Saving the mesh

save('sWT_1_Male_BuildedMesh.mat')