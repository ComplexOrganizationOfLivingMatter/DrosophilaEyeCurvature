function [uniformTrianglesCoordinates,uniformTriangles,P] = GetUniformMesh(baseIdenticalTriangles,heightIdenticalTriangles,coordinates,elements,k)

%% Inputs

% baseIdenticalTriangles: The initial guess of the base length of the
%                         uniform triangles.

% heightIdenticalTriangles: The initial guess of the height of the uniform
%                           triangles.

% coordinates: Matrix coordinates that represents the sWT mesh.

% elements: Connectivity matrix that represents coordinates.

% k: Vector whose components are the indexes of the points of coordinates
%    that belong to its boundary.

%% Outputs

% uniformTrianglesCoordinates: Coordinate tensor that represents the
%                              uniform mesh.

% uniformTriangles: Connectivity matrix that represents to the coordinate
%                   matrix of uniformTrianglesCoordinates

% P: Coordinate matrix that represented by uniformTriangles.

%% Code explanation

% This code will generate an uniform lattice in the x-y plane. We will set
% the boundary of sWT simulation and determine the points of the uniform
% lattice that are inside that boundary. If the number of triangles differs
% 5% from the number of triangles of sWT the algorithm will stop, in any
% other case, the uniform lattice would be changed (getting bigger or
% smaller) until the 5% criterio is satisfied.

%% Setting computational parameters

% We will apply this algorithm iteratively so in order to avoid
% computational cost (conditions as while) we will set a number of
% iterations high enough to be sure that the 5% condition is satisfied.

iterationNumber = 1000; % it will be far fewer itarions.

% Also, the dimensions for the uniform lattice (number of rows and columns)
% are set manually to, again, have a big enough region to avoid the sWT
% boundary to have empty spaces inside it. The only restriction is that
% they can be and odd number.

numberTrianglesRectangleBase = 100;
numberTrianglesRectangleHeight = 100;

%% Algorithm to obtain the uniform mesh

% Now, we will apply the iterative algorithm

for l = 1:iterationNumber

    % We start generation the base and the height of the uniform lattice:

    basePoints = linspace(0,(numberTrianglesRectangleBase + 1)*baseIdenticalTriangles,numberTrianglesRectangleBase);    
    basePointsDisplaced = basePoints + baseIdenticalTriangles/2;    
    heightPoints = linspace(0,(numberTrianglesRectangleHeight + 1)*heightIdenticalTriangles,numberTrianglesRectangleHeight);
    
    % It's important to remark that we add 1 because we need at least 2
    % points to define a triangle in both dimensions.
    
    % Now, we obtain the coordinate matrix of the uniform lattice:

    uniformLatticeCoordinates = GetUniformLatticeCoordinates(basePoints,basePointsDisplaced,heightPoints);
    
    % The segmentation data could not match with the region in which we 
    % have created the uniform lattice so we will match the barycenter of
    % those regions:
    
    sWTPoligon = polyshape(coordinates(k,1),coordinates(k,2));    
    [xsWTBarycenter,ysWTBarycenter] = centroid(sWTPoligon);    
    uniformLatticeBoundaryIndex = boundary(uniformLatticeCoordinates(:,1),uniformLatticeCoordinates(:,2));
    uniformLatticePoligon = polyshape(uniformLatticeCoordinates(uniformLatticeBoundaryIndex,1),uniformLatticeCoordinates(uniformLatticeBoundaryIndex,2));
    [xUniformLatticeBarycenter,yUniformLatticeBarycenter] = centroid(uniformLatticePoligon);
    
    % So the displacement vector is and the new sWT displaced boundary are:
    
    displacementVector = [xUniformLatticeBarycenter-xsWTBarycenter,yUniformLatticeBarycenter-ysWTBarycenter];
    newsWTBoundary = coordinates(k,:) + displacementVector;

    % Now we decide which points are inside the new sWT boundary:
    
    insidePoints = inpolygon(uniformLatticeCoordinates(:,1),uniformLatticeCoordinates(:,2),newsWTBoundary(:,1),newsWTBoundary(:,2));
    insideUniformLatticePoints = uniformLatticeCoordinates(insidePoints,:);
    
    % Now we triangulate all the space and impose the conditions:
    
    S = alphaShape(insideUniformLatticePoints);    
    S.Alpha = 1.1*S.Alpha;    
    [Tri,P] = alphaTriangulation(S);

    % Finally we impose the conditions to stop or to iteratively apply the
    % algorithm:
    
    if size(Tri,1) <= 1.01*2*size(elements,1) && size(Tri,1) >= 0.99*2*size(elements,1)
       
        break

    else

        if size(Tri,1) > 2*size(elements,1)
            
            baseIdenticalTriangles = 1.05*baseIdenticalTriangles;

        else

            baseIdenticalTriangles = 0.95*baseIdenticalTriangles;

        end

    end

end

%% Transform the uniform mesh into the set of its validTriangles

% Once we have the coordinate matrix and the connectivity matrix of the
% uniform lattice we will transform them to the tensor form for coherence:

uniformCoordinatesTensor = getCoordinatesTensor(P,Tri);

% Now, we can chose the corresponding set of validTriangles:

[uniformTrianglesCoordinates,uniformTriangles] = GetUniformTriangulation(uniformCoordinatesTensor,P,Tri);

