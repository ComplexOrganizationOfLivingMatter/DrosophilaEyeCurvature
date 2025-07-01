function [ConnectivityMatrixSmooth,CoordinatesSmooth] = LaplacianEnergySmoothingPupa(ConnectivityMatrix,Coordinates,lambda,iterations)

% This code will smooth the triangulation using the Laplacian Energy
% Method. The we will triangulate again the points:

Coordinates_Eye_Surface_Laplacian_Smoothed = LaplacianSmoothing(Coordinates,ConnectivityMatrix,lambda,iterations);

% Finally we triangulate the points:

tri2D = delaunay(Coordinates_Eye_Surface_Laplacian_Smoothed(:,1),Coordinates_Eye_Surface_Laplacian_Smoothed(:,2));
dt = triangulation(tri2D,Coordinates_Eye_Surface_Laplacian_Smoothed(:,1),Coordinates_Eye_Surface_Laplacian_Smoothed(:,2),Coordinates_Eye_Surface_Laplacian_Smoothed(:,3));
CoordinatesSmooth = dt.Points;
ConnectivityMatrixSmooth = dt.ConnectivityList;