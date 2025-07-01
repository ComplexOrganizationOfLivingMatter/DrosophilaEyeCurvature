function trianglesSetCoordinatesTensor = getCoordinatesTensor(P,trianglesSet)

%% Input

% P: Coordinate matrix of Nx2, being N the number of points.

% trianglesSet: Connectivity matrix that represent the set of triangles
%               whose coordinates in tensor form we want to obtain. It's a
%               Mx3 matrix, being M the number of triangles.

%% Output

% trianglesSetCoordinatesTensor: Coordinate tensor of Mx3x2. Each row 
%                                represents the 3 vertices of one triangle
%                                The x-y coordinates of each vertex are
%                                inside the 2 tensor dimension.

%% Code

trianglesSetCoordinatesTensor = zeros(size(trianglesSet,1),size(trianglesSet,2),size(P,2));

for i = 1:size(trianglesSetCoordinatesTensor,1)

    for j = 1:size(trianglesSetCoordinatesTensor,2)

        for k = 1:size(trianglesSetCoordinatesTensor,3)

            trianglesSetCoordinatesTensor(i,j,k)=P(trianglesSet(i,j),k);

        end

    end
    
end