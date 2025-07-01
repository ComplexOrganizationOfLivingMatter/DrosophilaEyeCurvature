function triangle_area = Triangle_Area_NonLinearTriangulation(elements,coordinates)

%% Inputs

% elements: It's a row (1 x 6) whose components are the rows of coordinates
%           that make up the triangle, which needs to have 6 points, 3 
%           vertices and 3 midpoints. The odd components are the vertices
%           and the even ones are the midpoints. elements needs to
%           refer to the points in clockwise or anti clockwise.

% coordinates: It's a matrix (n x 3, being n the number of points in the
%              entire triangulation) whose row i components are the
%              coordinates of the the point i in the triangulation.

%% Outputs

% triangle_area: The area of the triangle defined by the 6 points in
%                elements.

%% Code description

% We will calculate the area of the triangle defined by the 6 points (it's
% a non linear triangulation) that make up the triangle. To do that, we
% will split the triangle in 4 triangles, whose areas will be calculated
% and then summed up to obtain the total area.

% The only restriction in elements is the clockwise or anticlockwise, so
% this means that in one triangle the point the first component of elements
% refer to is the vertex of the top, while in other triangle it's one of
% the vertices of the base (isosceles triangle for example).

% This restriction is quite useful because it gives us the closest points.
% For example, the closest points to elements(1) are elements(2) and
% elements(6), the closest points to elements(3) are elementes(2) and
% elements(4)... With the vertices and the midpoints we will have 3
% triangles, the last fourth is the triangle maked up by the 3 midpoints.

%% Code

% We will calculate the area of the triangles one by one:

area_triangle_1 = Triangle3D_AreaCalculator(coordinates(elements([1,2,6]),1),coordinates(elements([1,2,6]),2),coordinates(elements([1,2,6]),3));

area_triangle_2 = Triangle3D_AreaCalculator(coordinates(elements([3,2,4]),1),coordinates(elements([3,2,4]),2),coordinates(elements([3,2,4]),3));

area_triangle_3 = Triangle3D_AreaCalculator(coordinates(elements([5,4,6]),1),coordinates(elements([5,4,6]),2),coordinates(elements([5,4,6]),3));

area_triangle_4 = Triangle3D_AreaCalculator(coordinates(elements([2,4,6]),1),coordinates(elements([2,4,6]),2),coordinates(elements([2,4,6]),3));

% Finally, we have that:

triangle_area = area_triangle_1 + area_triangle_2 + area_triangle_3 + area_triangle_4;










