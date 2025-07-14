function [alpha, beta] = unitVectorToAngle(u)

%% Input

% u: vector whose angles between the xy plane and the z-axis we want to
%    obtain.

%% Outputs

% alpha: Angle between u vector and xy plane.

% beta: Angle between u vector and z-axis.

%% Code explanation

% This code will obtain the angles of vector u between xy plane and the
% x-axis.

%% Rotational angle between the projected u on the xy plane and the x-axis

alpha = atan2(u(2), u(1)); 

%% Rotational angle between the u vector and the z-axis

beta = atan2(sqrt(u(1)^2 + u(2)^2), u(3));