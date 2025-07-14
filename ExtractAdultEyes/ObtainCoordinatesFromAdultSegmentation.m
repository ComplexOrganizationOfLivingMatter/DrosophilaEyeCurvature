function coordinates = ObtainCoordinatesFromAdultSegmentation(adultEyeSegmentationImage)

%% Input

% adultEyeSegmentationImage: 3D Image in .tif binarize resulting from the
%                            semantic segmentation of the adult eye.

%% Output

% coordinates: 3D coordinate matrix that represent the volume of 
%              adultEyeSegmentationImage in microns.

%% Code explanation 

% This code will read adultEyeSegmentationImage to extract the volumetric
% coordinates of the adult eye surface. Then, using the resolution in x,y,z
% the coordinates will be trasnformed to microns.

%% Code

% First we will read the adultEyeSegmentationImage and the x,y,z
% resolution:

[img,imgInfo] = readStackTif(adultEyeSegmentationImage);
xyResolution = 1/imgInfo(1).XResolution;
spacingInfo = strsplit(imgInfo(1).ImageDescription, 'spacing=');
spacingInfo = strsplit(spacingInfo{2}, '\n');
zResolution = str2double(spacingInfo{1});

% Now, we obtain the volumetric coordinates:

bwimg = double(bwperim(img));
[x,y,z] = ind2sub(size(bwimg),find(bwimg==1));

% Finally we transform then into microns:

x = x*xyResolution;
y = y*xyResolution;
z = z*zResolution;

coordinates = [x,y,z];

   
