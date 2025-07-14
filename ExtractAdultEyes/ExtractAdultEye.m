function [Coordinates_Eye_Surface,Connectivity_matrix] = ExtractAdultEye(adultEyeSegmentationImage,alphaFactor)

%% Inputs

% adultEyeSegmentationImage: 3D Image in .tif binarize resulting from the
%                            semantic segmentation of the adult eye.

% alphaFactor: Factor ranging from 0 (empty triangulation) to infinite
%              (convexHull). It determines the sensitivy at taking the eye
%              surface boundary. It's recommend a 1.8-2.2 value.

%% Outputs

% Coordinates_Eye_Surface: Coordinate matrix whith the coordinates of the
%                          3D adult eye.

% Connectivity_matrix: Connectivity matrix that represents the
%                      Coordinates_Eye_Surface coordinate matrix.

% REMARK: The nomenclature of the output (snake case) differs from the
% usual nomenclature (camel case). It's due to technical inheritance, so
% it's better not to change.

%% Code Explanation 

% This code will read the adultEyeSegmentationImage in order to extract
% it's x,y,z-resolution and the coordinates of the pixels that belong to
% the adult eye in microns. Then the adult eye will be oriented to lay in
% the x-y plane in such a way that the apex matches the point with highest
% z-coordinates. Afterwards, an algorithm to extract the coordinates of the
% eye surface will be applied.

% Remarking point: Extract the surface through image analysis using the
% function bwmorph(bw, 'thin', Inf) could lead to irregularities in the
% surface due to a miss match with the segmentation process.
 
%% Reading the image and extracting the coordinates in microns

% We will read the image and extract the coordinates in microns of the
% surface 

coordinates = ObtainCoordinatesFromAdultSegmentation(adultEyeSegmentationImage);

%% Orienting coordinates

% Now we will orient the coordinates matrix to match the apex with the
% point point of highest z-value. We first rotate the coordinates in such a
% way the volumetric coordinates are align its major azis to z-axis.

coordinates = rotateCoordinatesAlongZ(coordinates,'x');

% To choose the adult eye surface (the exterior surface of the eye) we need
% the z-axis to be "normal" to the surface. This means that we need to
% rotate the eye around x-axis pi/2. To avoid more corrections, we will set
% the point with the minimun z-coordinate to the origin of coordinates.
% This way, we will rotate the eye but not this point, so, after the
% rotation, the eye will be set in z = 0 as the minimum z-coordinate: 

zMinIndex = find(coordinates(:,3) == min(coordinates(:,3)));
displacementVector = -coordinates(zMinIndex(1),:);
coordinates = coordinates + displacementVector;

% Now, we rotate the eye around the x-axis pi/2 (counter-clockwise). To do
% that, first, we define the coordinates as a pointcloud object:

eyeCloud = pointCloud(coordinates);
theta= -90; %counter-clockwise
rot = rotx(theta);
trans = [0,0,0];
tform = rigid3d(rot,trans);
eyeCloudRotated = pctransform(eyeCloud,tform);
coordinates = eyeCloudRotated.Location;

%% Taking the boundary: Interpolating and triangulating the coordinates.

% Since we need to measure to adimensional gaussian curvature over the
% adult eye surface we need to capture the boundary as good as posible so
% we need to interpolate the volume.

coordinates = InterpolateCoordinates(coordinates);

% Now, we use alphaShape function to triangulate the coordinates. It's
% importante to remark that we are interested in the boundary so using a
% delaunay triangulation could lead to artifacts in the boundary since
% delaunay triangulation is a convex triangulation.

S = alphaShape(coordinates(:,1),coordinates(:,2),coordinates(:,3));
S.Alpha = alphaFactor*S.Alpha;
coordinates = S.Points;

% Now we obtain the boundary of coordinates in the xy plane:

boundaryPoints = ObtainBoundaryPoints(coordinates);

%% Applying the algorithm: Extracting the adult eye surface

% Now the algorithm is clear, we start at the point of maximum z and we
% keep al the triangles connected to him. We repeat that with the new
% points (triangles connected to the points in the previous iteration)
% until we reach the triangles whose vertices are at least, one of the
% points that appear un boundary_rows:

coordinatesEyeSurface = ObtainAdultEyeFromVolume(coordinates,boundaryPoints,S);

%% Having a good triangulation and getting the outputs

% First, we inteprolate again the coordinates to have a "continuous" mesh:

coordinatesEyeSurfaceInterpolated = InterpolateCoordinates(coordinatesEyeSurface);

% If we use alphashape to triangulate the surface it can have holes and
% other kind of computational artifacts. So, it's better to use a delaunay
% triangulation in 2D and then use the triangulation function to obtain the
% 3D triangulation.

tri = delaunay(coordinatesEyeSurfaceInterpolated(:,1),coordinatesEyeSurfaceInterpolated(:,2));
tri3D = triangulation(tri,[coordinatesEyeSurfaceInterpolated(:,1),coordinatesEyeSurfaceInterpolated(:,2),coordinatesEyeSurfaceInterpolated(:,3)]);

% Finally, we obtain the surface points and the triangulation:

Coordinates_Eye_Surface = tri3D.Points;
Connectivity_matrix = tri3D.ConnectivityList;
