function [validTriangles,inverseTriangles,P,validTrianglesCoordinates,inverseTrianglesCoordinates,conversionFactor] = ExtractTriangulationFromPupa(fileName,path2save,initialTriangle,alphaFactor)

%% Inputs

% fileName: Name of the file which contains the labels from the
%           segmentation. It must be a .mat file.

% path2save: Path to save the results in .mat.

% initialTriangle: Index of the inital triangle to start the algorithm.

% basalImage: Basal image of the corresponding pupa. It must be a 2D .tif
%             with the same size as the labels of fileName.

% alphaFactor: Number to multiply S.Alpha in the triangulation to manually
%              correct the mesh.

%% Outputs

% validTriangles: Connectivity matrix that represents the triangles that
%                 correspond to the set of valid triangles.

% inverseTriangles: Connectivity matrix that represents the triangles that
%                   correspond to the set of invalid triangles.

% P: Coordinate matrix corresponding to the set of valid triangles and the
%    inverse triangles.

% validTrianglesCoordinates: Coordinate tensor wich contains the
%                            coordinates of the validTriangles. It's a
%                            tensor because it makes easier the
%                            postprocessing.

% inverseTrianglesCoordinates: Coordinate tensor wich contains the
%                              coordinates of the inverseTriangles. It's a
%                              tensor because it makes easier the
%                              postprocessing.

% conversionFactor: It's the conversion factor from pixels² to microns².
%                   It's the square power of the spatial revolution of 
%                   basalImage.

%% Code Explanation

% This code will load the file under the name fileName, which contains the
% labels from the segmentation of the corresponding pupa (with the same
% fileName). Then, it will obtain the total mesh of the labels using an
% alpha triangulation adn taking the initialTriangle it will obtain the
% validTriangles and inverseTriangles sets of triangles. Finally, it will
% plot the set of validTriangles over the basalImage.

%% Load labels, transform to double and extract coordinates

% The labels are expected to be logical:

load(filaName,'labels');
labels  = double(sum(labels,3));

% In order to use regionprops3 we need to label the diferrent elements,
% because we need the different centroids. So, we will use bwlabel(img) which
% assigns labels to different elements. For example:

% img = [1,0,0;0,1,0,0,0,1]; it has 3 different elements.

% If the connection is 0 (by default is 8 and it can be change to 4 using
% that bwlabel(img,connection) is a property of bwlabel.

% bwlabel(img) = [1,0,0;,0,2,0;0,0,3];

differentLabels = bwlabel(labels);

% Now we obtain all its properties using regionprops3, specially, the
% coordinates of differentLabels.

geometricProperties = regionprops3(differentLabels);
centroidTable = geometricProperties.Centroid;

%% Building the mesh and the  growth algorithm

% As our tissue is not convex we will use alphaShape(P) where P is a matrix
% with the points in the triangulation. 

S = alphaShape(centroidTable);

% In order to fill all the space we need to modify the S.Alpha factor (the
% greater it is, filler is the space).

S.Alpha = alphaFactor*S.Alpha;

% We need the vertices that form every triangle and its coordinates. The
% vertices are in Tri matrix (row i es triangle i adn the number are the ID
% of vertex) and the coordinates in P (row i are the coordinates of vertex
% i)

[Tri,P] = alphaTriangulation(S);

% We apply the region-growing algorithm: 

[validTriangles,inverseTriangles] = regionGrowingAlgorithm(Tri,initialTriangle);

%% Obtaining validTrianglesCoordinates and inverseTrianglesCoordinates

validTrianglesCoordinates = getCoordinatesTensor(P,validTriangles);
inverseTrianglesCoordinates = getCoordinatesTensor(P,inverseTriangles);

%% Obtaining area from validTriangles previous to its representation

% We need to obtain the conversion factor from pixels to microns using the
% image info:

imgInfo = imfinfo(basalImage);
conversionFactor = 1/imgInfo.XResolution^2;

% Now, we obtain the area distribution of validTriangles:

validTrianglesAreaDistribution = getAreaDistributionFromTriangulation(validTrianglesCoordinates,conversionFactor);

%% Initialize the figure with the basalImage

% To plot the set of validTriangles over the basalImage we need to load the
% basalImage.

photo = imread(basalImage);

% Recommended: Use the photo in uint16 format instead of uint8 and adjust
% the intensity values:

photo = imadjust(uint16(photo));

% The figure is:

figure
imshow(photo)
caxis('auto');
colormap
hold on

%% Setting the parameters and colormap of the figure

% The minimum and maximum area to normalize all pupas are:

AreaMin = 53.0170;  % Don't change
AreaMax = 249.8067; % Don't change

% The colorblind scale is (don't change):

colors = [ ...
    0.00, 1.00, 0.00;  % Green
    1.00, 0.00, 1.00;  % Magenta
    0.00, 1.00, 1.00;  % Cyan
    1.00, 1.00, 1.00]; % White

% Number of interpolation steps (smoother = more steps)
nSteps = 256;

% Interpolate the colors to create a smooth gradient
colors_area = interp1(linspace(0, 1, size(colors, 1)), colors, linspace(0, 1, nSteps));
AreaIntervals = linspace(AreaMin,AreaMax,size(colors_area,1));

% Finally we change the colormap of the figure:

d = double(photo);
grayscaleIM(:,:,3) = (d-min(d(:)))/(max(d(:))-min(d(:)));
grayscaleIM(:,:,2) = (d-min(d(:)))/(max(d(:))-min(d(:)));
grayscaleIM(:,:,1) = (d-min(d(:)))/(max(d(:))-min(d(:)));
colormap(colors_area)
imagesc(grayscaleIM)

%% Plotting the mesh over the basal image

for i = 1:size(validTriangles,1)

        pgon = polyshape(validTrianglesCoordinates(i,:,1),validTrianglesCoordinates(i,:,2));

    for k = 1:length(AreaIntervals)-1

        if AreaIntervals(k) <= validTrianglesAreaDistribution(i) && validTrianglesAreaDistribution(i) <= AreaIntervals(k+1)

            plot(pgon,'FaceColor',colors_area(k,:),'FaceAlpha',1)
            hold on

        end

    end

end

caxis([AreaMin,AreaMax])
axis off

%% Saving the results

name2saveVariables = strrep(fileName,'_Labels.mat','_MeshVariables.mat');
name2saveFigure = strrep(fileName,'_Labels.mat','_BasalMesh.fig');

save(strcat(path2save,name2saveVariables),'validTriangles','inverseTriangles','P','validTrianglesCoordinates','inverseTrianglesCoordinates')
savefig(strcat(path2save,name2saveFigure))
