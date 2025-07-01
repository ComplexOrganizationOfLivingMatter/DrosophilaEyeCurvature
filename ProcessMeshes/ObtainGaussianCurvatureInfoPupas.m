function ObtainGaussianCurvatureInfoPupas(pupasPath,outputPathCurvature)

% This code will obtain the triangulations from the stl files in pupasPath
% and store all the curvature information in outputPathCurvature and the
% tables in outputPathTable.

pupasDir = dir(strcat(pupasPath,'*','.stl'));
n = 10;

for i = 1:size(pupasDir,1)

    % First, we read the pupas and obtain the variables:

    pupaName = pupasDir(i).name;
    triangulation = stlread(strcat(pupasPath,pupaName));
    Coordinates_Eye_Surface_Regularized = triangulation.Points;
    Connectivity_matrix_Regularized = triangulation.ConnectivityList;

    % Now, we obtain all the curvature info:

    [Coordinates_Eye_Surface_Regularized_Interior,K_Interior] = compute_gaussian_curvature(Coordinates_Eye_Surface_Regularized,Connectivity_matrix_Regularized);
    [GaussianCurvaturePercentiles,Coordinates_Eye_Surface_Regularized_Interior,points_in_percentil,K_Interior] = GetGaussianCurvaturePercentiles(Coordinates_Eye_Surface_Regularized_Interior,K_Interior,n);
    pupaName2save = strrep(pupaName,'.stl',strcat('_GCD_n_',num2str(n),'.mat'));
    save(strcat(outputPathCurvature,pupaName2save),'Coordinates_Eye_Surface_Regularized','Connectivity_matrix_Regularized','Coordinates_Eye_Surface_Regularized_Interior','K_Interior',"GaussianCurvaturePercentiles","points_in_percentil")

end

