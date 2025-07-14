function u = pcaEig(pc, magnitude)

%% Inputs

% pc: Coordinate matrix (called pc for point cloud but not being a point
%     cloud object) whose eigenvector of a certain magnitude we want to
%     obtain.

% magnitude: String that defines the eigenvector you want to obtain. 'max'
%            is for highest eigenvalue, 'middle' for the second and 'min' 
%            for the lowest ones.

%% Output

% u: eigenvector of pc defined by magnitude

%% Code explanation

% This code will obtain the eigenvector of pc using the covariance matrix
% and normal vectors.

%% Obtain the covariance matrix

covariance = cov([pc(:, 1) pc(:, 2) pc(:, 3)]);

%% Calculate the eigenvectors and obtain the normal vector

[V, D] = eig(covariance);

diagonalEigenvalues = diag(D);

%% Output the normal vector 

% Sort the eigenvectors based on size of eigenvalues 

[~, I] = sort(diagonalEigenvalues);
V = V(:, I);
    
switch magnitude
    case 'max'
        % Choose the eigenvector of the highest eigenvalue
        u = V(:, 3); 
    case 'middle'
        % Choose the eigenvector of the middle eigenvalue
        u = V(:, 2); 
    case 'min'
        % Choose the eigenvector of the lowest eigenvalue
        u = V(:, 1); 
end