function pointloadFixed = FixPressurePoinload(pointload,normals,originalLoad)

%% Inputs

% pointload: The first column is the ID of the point, the second column is
%            the degree of freedom of that point and the third column is
%            the applied load to that point in that degree of freedom.

% normals: The row i is the normal vector to the point pointload(i,1).

% originalLoad: Scalar. Is the value of the pointload applied in the
%               first iteration.

%% Outputs

% pointloadFixed: The pointload matrix now with the first 3 degree of
%                 freedom of each node being applied a pressure given by
%                 normals and originalLoad.
%% Code explanation

% The matrix pointload just have the third degree of freedom (z direction)
% being applied a load. As we want to applied a pressure, which is always
% normal to the point, we need to descompose the original load in the x and
% y direction. These new load must be added to pointload matrix in
% pointloadFixed.

%% Code

% First, we obtain the ID of the points in poinload and we obtain 3 sorted
% copies of them because the pressure will affect to the 3 first degree of
% freedom of every point:

pointsID = repmat(pointload(:,1),3,1);
pointsID = sort(pointsID);

% Now, for every point, we want the degrees of freedom, so we need
% size(pointload,1) copies of the column [1;2;3]:

pointsFreedomDegree = repmat([1;2;3],size(pointload,1),1);

% Now, we multiply the normals by originalLoad to obtain:
normalsIDbyID = originalLoad*reshape(normals.',1,[])';

% Finally, we obtain the pointloadFixed:

pointloadFixed = [pointsID,pointsFreedomDegree,normalsIDbyID];




