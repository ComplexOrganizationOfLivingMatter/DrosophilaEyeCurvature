%% Main Process Pupas

pupasPath = '/media/juan/ce04cf3c-0227-4dac-9535-9f014799bf45/juan/Projects/Auxetics/Auxetic_Eyes/Eyes_Expansion_Codes/Analysis/Pressure/D_Mauritiana_Triangulation_RedDots/mauritaniaTam16 male42h25 retina12 actin/Pupa_Iterations/';
outputPath = '/media/juan/ce04cf3c-0227-4dac-9535-9f014799bf45/juan/Projects/Auxetics/Auxetic_Eyes/Eyes_Expansion_Codes/Analysis/Pressure/D_Mauritiana_Triangulation_RedDots/mauritaniaTam16 male42h25 retina12 actin/Pupa_Iterations_Inter_DS_LES_stl/';

iterations = 75;
lambda = 0.1;

ProcessPupas_Inter_DS_LES_STL(pupasPath,outputPath,iterations,lambda)

% iterations = 75 for all pupas except IT; 750 for IT pupas.
% lamba = 0.1 for all pupas except IT; lambda = 0.01 for IT pupas.
