function Triangle_Area = Triangle3D_AreaCalculator(x,y,z)

%% Inputs

% x: The x-coordinates of the vertices. It's a 3x1 column vector.

% y: The y-coordinates of the vertices. It's a 3x1 column vector.

% z: The z-coordinates of the vertices. It's a 3x1 column vector.

%% Outputs

% Triangle_Area: The area of the triangle

%% Code
x = x';
y = y';
z = z';
ons = [1 1 1];
Triangle_Area = 0.5*sqrt(det([x;y;ons])^2 + det([y;z;ons])^2 + det([z;x;ons])^2);