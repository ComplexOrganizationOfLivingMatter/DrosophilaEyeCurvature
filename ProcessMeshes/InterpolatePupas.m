function [coordinatesMatrix] = InterpolatePupas(coordinates)

% This code will interpolate the eyes to obtain a dense point cloud of
% them:

x = coordinates(:,1);
y = coordinates(:,2);
z = coordinates(:,3);
[xq,yq] = meshgrid(min(x):max(x), min(y):max(y));
zq = griddata(x,y,z,xq,yq);
x = xq(:);
y = yq(:);
z = zq(:);
coordinatesMatrix = [x,y,z];
coordinatesMatrix(isnan(z),:) = [];

