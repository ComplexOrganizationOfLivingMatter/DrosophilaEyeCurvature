function interpolatedCoordinates = InterpolateCoordinates(coordinates)

%% Input

% coordinates: 3D coordinate matrix that represent the volume of 
%              adultEyeSegmentationImage in microns.

%% Output

% interpolatedCoordinates: 3D coordinate matrix that represent the volume 
%                          of adultEyeSegmentationImage in microns
%                          interpolated.

%% Code explanation

% This code will interpolate the volume of the 3D coordinate matrix
% coordinates.

%% Code

% First, we define the x, y and z coordinates:

x = coordinates(:,1);
y = coordinates(:,2);
z = coordinates(:,3);

% Now, using meshgrid we interpolate over the plane xy:

[xq,yq] = meshgrid(min(x):max(x), min(y):max(y));

% Once we have the interpolation over xy plane, we interpolate using
% griddata to obtain the interpolated z coordinates.

zq = griddata(x,y,z,xq,yq);
x = xq(:);
y = yq(:);
z = zq(:);

% Finally, we obtain the inteporlateCoordinates matrix. As coordinates is a
% dense volumen the interpolation may fail ar some points, givin a nan in
% the z coordinate. These nans are removed:

interpolatedCoordinates = [x,y,z];
interpolatedCoordinates(isnan(z),:) = [];