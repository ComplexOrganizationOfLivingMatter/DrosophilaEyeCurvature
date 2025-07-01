function normals = findNormals(coordinates)

%% Input

% coordinates: Coordinates matrix as n x 3, being n the number of points.
%              The points must belong to a 3d surface.

%% Output

% normals: Normal vector to every points specified as a n x 3 matrix, being
%          n the number of points. The row i is the unit normal vector to 
%          the row i of coordinates

%% Code explanation

% This code will obtain the unit normal vectors to the points in
% coordinates. The unit normal vector can point inwards or outwards the
% surface, and we want the vector to point outwards. To be sure of that, we
% will obtain the centroid in x-y plane, obtain the position vector of the
% points in coordinates from that centroid and will check if the normal
% vector and the position vector have the same direction.

%% Code

% First we obtain the centroid in the plane:

planeCentroid = [mean(coordinates(:,1)),mean(coordinates(:,2)),min(coordinates(:,3))];

% Now, we transform the coorduinates in a pointCloud and obtain the
% normals:

ptCloud = pointCloud(coordinates);
normals = pcnormals(ptCloud);
points = ptCloud.Location;

% Now, we do the loop over every point and fix the normals:

for i = 1:size(points,1)
    
    % We compute the dot product of the normal and the vector from the 
    % planeCentroid:

    pointVector = points(i,:) - planeCentroid;
    dot_product = dot(normals(i,:), pointVector);
    
%     If the dot product is negative, flip the normal:

    if dot_product < 0
        normals(i,:) = -normals(i,:);
    end

    if normals(i,3) < 0 

        normals(i,:) = -normals(i,:);

    end
    
end

