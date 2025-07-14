function pc2 = rotatePC(pc, Ry, Rz)

%% Inputs

% pc: Coordinate matrix (called pc for point cloud but not being a point
%     cloud object) that we want to rotate.

% Ry: Rotational matrix around y-axis.

% Rz: Rotational matrix around z-axis.

%% Output

% pc2: Coordinate matrix (called pc2 for point cloud but not being a point
%     cloud object) that is the result of rotatin pc using Ry and Rz.

%% Code explanation

% This code will rotate the pc coordinate matrix using the rotational
% matrices Ry and Rz.

%% Code

% Convert the point cloud to 3 * N format and rotate around the z-axis to
% align point cloud along x-axis:

matrix = pc';
matrix2 = Rz*matrix;
matrix2 = Ry*matrix2;

% Finally, the rotated point cloud is:

pc2 = matrix2';