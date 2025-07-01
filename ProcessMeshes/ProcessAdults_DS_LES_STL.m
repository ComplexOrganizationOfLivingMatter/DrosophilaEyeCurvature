function ProcessAdults_DS_LES_STL(adultsPath,outputPath,iterations)

%% Pipeline to obtain pupas in stl format to meshlab

% This code will downsample the adult eyes pointCloud, then we will apply
% the Laplacian Energy Smoothing and write it in a stl file wich will be 
% then processed.

%% Code

adultsDir = dir(strcat(adultsPath,'*','.mat'));

for i = 1:size(adultsDir,1)

    adultName = adultsDir(i).name;
    load(strcat(adultsPath,adultName),'Coordinates_Eye_Surface','Coordinates_Eye_Surface')

    % Now we downsample the pointCloud:

    [Connectivity_matrix_DS, Coordinates_DS] = DownSamplePupas(Coordinates_Eye_Surface);

    % Now we smooth the eyes using the Laplacian Energy Smoothing method:

    lambda = 0.1;

    [Connectivity_matrix_DS_LES,Coordinates_DS_LES] = LaplacianEnergySmoothingPupa(Connectivity_matrix_DS,Coordinates_DS,lambda,iterations);
    
    % Finally, sa save the results in a stl file:

    TR = triangulation(Connectivity_matrix_DS_LES,Coordinates_DS_LES(:,1),Coordinates_DS_LES(:,2),Coordinates_DS_LES(:,3));
    adultName2save = strrep(adultName,'.mat','_DS_LES.stl');
    stlwrite(TR,strcat(outputPath,adultName2save));

end



