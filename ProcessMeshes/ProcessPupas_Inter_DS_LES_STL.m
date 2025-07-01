function ProcessPupas_Inter_DS_LES_STL(pupasPath,outputPath,iterations,lambda)

%% Pipeline to obtain pupas in stl format to meshlab

% This code will interpolate the pupas, them will do the downsample of the
% pointCloud, then we will apply the Laplacian Energy Smoothing and write
% it in a stl file wich will be then processed.

%% Code

pupasDir = dir(strcat(pupasPath,'*','.mat'));

for i = 1:size(pupasDir,1)

    pupaName = pupasDir(i).name;
    load(strcat(pupasPath,pupaName),'coordinates','elements')
    goodelements = [elements(:,1),elements(:,3),elements(:,5)];
    goodelements = unique(goodelements);
    coordinates = coordinates(goodelements,:);

    % Now we interpolate the eye:

    [Coordinates] = InterpolatePupas(coordinates);

    % Now we downsample the pointCloud:

    [Connectivity_matrix_DS, Coordinates_DS] = DownSamplePupas(Coordinates);

%     tr = triangulation(Connectivity_matrix_DS,Coordinates_DS(:,1),Coordinates_DS(:,2),Coordinates_DS(:,3));
%     pathAuxiliar = '/media/juan/ce04cf3c-0227-4dac-9535-9f014799bf45/juan/Projects/Auxetics/Auxetic_Eyes/Reading_3D_Eyes/Calibration/Results/SmothedEyes/Pupas_Mat_DS/';
%     save(strcat(pathAuxiliar,pupaName),'tr');

    % Now we smooth the eyes using the Laplacian Energy Smoothing method:

    [Connectivity_matrix_DS_LES,Coordinates_DS_LES] = LaplacianEnergySmoothingPupa(Connectivity_matrix_DS,Coordinates_DS,lambda,iterations);
    
    % Finally, sa save the results in a stl file:

    TR = triangulation(Connectivity_matrix_DS_LES,Coordinates_DS_LES(:,1),Coordinates_DS_LES(:,2),Coordinates_DS_LES(:,3));
    pupaName2save = strrep(pupaName,'.mat','_Inter_DS_LES.stl');
    stlwrite(TR,strcat(outputPath,pupaName2save));

end



