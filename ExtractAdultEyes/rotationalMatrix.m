function [Rx, Ry, Rz] = rotationalMatrix(alpha, beta)

%% Inputs

% aplha: Angle to rotate the coordinates around x-axis.

% beta: Angle to rotate the coordinates around y-axis.

%% Ouputs

% Rx: Rotational matrix around x-axis.

% Ry: Rotational matrix around y-axis.

% Rz: Rotational matrix around z-axis.

%% Code explanation

% This code will build the rotational matrices around x, y and z axes.

%% Code
    
Rx = [1 0 0; 0 cos(beta) -sin(beta); 0 sin(beta) cos(beta)];
Ry = [cos(beta) 0 sin(beta); 0 1 0; -sin(beta) 0 cos(beta)];
Rz = [cos(alpha) -sin(alpha) 0; sin(alpha) cos(alpha) 0; 0 0 1];