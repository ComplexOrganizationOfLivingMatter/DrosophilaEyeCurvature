%% Code Explanation

% This code will build the the rubber-like mesh version of sWT simulations.
% As we did for sWT simulations, this mesh will accomplish the requirements
% of the computation model, e.g, being  quadratic triangulation, define the
% Dirichlet Boundary Conditions and the Nodal points.

%%  Material Properties

% We define the porperties of the isometric material we will simulate and
% the pressure apllied to it;

young = 6.825e9;
poiss = 3e-1;
denss = 0;
thick = 0.125;
p = 1e5;

%% Building the rubber-like model from sWT simulation

% To obtain the rubber-like model we just need to load the boundary of sWT:

load("sWT_1_Male_BuiltMesh.mat",'coordinates','k');

% Now, we can obtain the the rubber-like model with coordinates and
% elements satisfying the requirements of the computational model:

[coordinates,elements] = GetRubberMesh(coordinates,k);

%% Dirichlet Boundary Conditions and nodal points

% As we are working  with a continuous material it doesn't make sense to
% choose the vertices of the triangles as the boundary. We will choose also
% de midpoints and we won't have any problem since the triangles are small
% enough to not be deformated by the tension due to the midpoints. So, we
% have that:

k = boundary(coordinates(:,1:2),0.99);

% Once the boundary has been identified, the Dirichlet Boundary Conditions
% are defined in fixnodes and the nodal points in pointload:

fixnodes = [k,ones(size(k,1),1),zeros(size(k,1),1);k,2*ones(size(k,1),1),zeros(size(k,1),1);k,3*ones(size(k,1),1),zeros(size(k,1),1)];
pointloadnodes = setdiff(unique(elements(:,[1,3,5])),k);
pointload = [pointloadnodes,3*ones(size(pointloadnodes,1),1),p*ones(size(pointloadnodes,1),1)];

% Side loads

uniload = sparse(size(elements,1),1);

%% Saving the mesh

save('sRub_1_Male_BuiltMesh.mat')