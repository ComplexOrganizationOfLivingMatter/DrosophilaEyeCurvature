%% Main Gaussian Curvature Info & minimum difference

pupasPath = '/media/juan/ce04cf3c-0227-4dac-9535-9f014799bf45/juan/Projects/Auxetics/Auxetic_Eyes/Eyes_Expansion_Codes/Analysis/Pressure/D_Mauritiana_Triangulation_RedDots/mauritaniaTam16 male42h25 retina12 actin/Pupa_Iterations_Inter_DS_LES_stl_IER/';
outputPathCurvature = '/media/juan/ce04cf3c-0227-4dac-9535-9f014799bf45/juan/Projects/Auxetics/Auxetic_Eyes/Eyes_Expansion_Codes/Analysis/Pressure/D_Mauritiana_Triangulation_RedDots/mauritaniaTam16 male42h25 retina12 actin/Pupa_Gaussian_Distribution/';
adultPath = '/media/juan/ce04cf3c-0227-4dac-9535-9f014799bf45/juan/Projects/Auxetics/Auxetic_Eyes/Reading_3D_Eyes/Calibration/Results/SmothedEyes/Gaussian_Curvature_Distribution/Definitive_Mauritiana_Good/';
% adultPathHausdorff = '/media/juan/ce04cf3c-0227-4dac-9535-9f014799bf45/juan/Projects/Auxetics/Auxetic_Eyes/Reading_3D_Eyes/Calibration/Results/SmothedEyes/Gaussian_Curvature_Distribution/Definitive_WT_Good_Oriented/';
outputPath = '/media/juan/ce04cf3c-0227-4dac-9535-9f014799bf45/juan/Projects/Auxetics/Auxetic_Eyes/Eyes_Expansion_Codes/Analysis/Pressure/D_Mauritiana_Triangulation_RedDots/mauritaniaTam16 male42h25 retina12 actin/Min_Difference_Iteration/';

% Obtain the curvature info and save in .mat:

ObtainGaussianCurvatureInfoPupas(pupasPath,outputPathCurvature);

% Obtain comparison table using the curvature info and Hausdorff:

[nGauss,minDifferenceGauss,minDifferenceIterationsGauss,pupaMeanDifferenceIterationsGauss] = ObtainMinimumDifference(outputPathCurvature,adultPath,outputPath);
% [nHausdorff,minDifferenceHausdorff,minDifferenceIterationsHausdorff,pupaMeanDifferenceIterationsHausdorff] = ObtainMinimumDifferenceHausdorff(outputPathCurvature,adultPathHausdorff,outputPath);

% We finally save the minDifference and n. First for the Gauss comparison:

pupasName = strsplit(outputPathCurvature,'/');
figure
plot(minDifferenceIterationsGauss)
xlabel('Iterations')
ylabel('Gauss Metric')
figTitle = strcat(pupasName{end-2},' Gauss Metric along iterations');
title(figTitle)
savefig(strcat(outputPath,pupasName{end-2},'_GaussMetric_Iterations.fig'));
save(strcat(outputPath,pupasName{end-2}),'minDifferenceIterationsGauss','minDifferenceGauss','nGauss','pupaMeanDifferenceIterationsGauss')

% Then for the Hausdorff comparison:

% figure
% plot(minDifferenceIterationsHausdorff)
% xlabel('Iterations')
% ylabel('Hausdorff Distance (\mum)')
% figTitle = strcat(pupasName{end-2},' Hausdorff Distance along iterations');
% title(figTitle)
% savefig(strcat(outputPath,pupasName{end-2},'_HausdorffDistance_Iterations.fig'));
% % save(strcat(outputPath,pupasName{end-2}),'minDifferenceIterationsHausdorff','minDifferenceHausdorff','nHausdorff','pupaMeanDifferenceIterationsHausdorff')
% save(strcat(outputPath,pupasName{end-2}),'minDifferenceIterationsGauss','minDifferenceGauss','nGauss','pupaMeanDifferenceIterationsGauss','minDifferenceIterationsHausdorff','minDifferenceHausdorff','nHausdorff','pupaMeanDifferenceIterationsHausdorff')

