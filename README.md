# DrosophilaEyeCurvature

## Introduction

This repository contains the code for the analysis of pupae and adult segmentation, build the meshes for the computational model, the computational model and to process the meshes. Follow the steps below to run or reproduce our analysis.

---

## Steps to Reproduce the Analysis

### 1. Extract the mesh from the pupae segmentation

Code in: TriangularMeshGenerationFromSegmentation/

Code for extracting the triangular mesh in the basal surface of Drosophila pupae. 

The segmentation must be load in .mat and the results will be saved in .mat as well.

Run [ExtractTriangulationFromPupa.m](TriangularMeshGenerationFromSegmentation/ExtractTriangulationFromPupa.m)
 to obtain the set of validTriangles, the set of inverseTriangles and their coordinates.

---

### 2. Build the mesh for the Computational Model

Code in: BuildMeshToComputaionalModel/

Code for obtaining the quadratic triangular mesh from the pupae basal triangular mesh.

Once the basal triangulation is obtained, run the different .m files in this folder to get them in the format that the computational model needs.

Run [Geometry_sWT_1_Male.m](BuildMeshToComputaionalModel/Geometry_sWT_1_Male.m) to obtain the sWT computational model information. Run [Geometry_sIT_1_Male.m](BuildMeshToComputaionalModel/Geometry_sIT_1_Male.m) and [Geometry_sRub_1_Male.m](BuildMeshToComputaionalModel/Geometry_sRub_1_Male.m) for its computational versions.

The data must be saved as .mat in .m file such as [sWT_1_Male_ComputationalModelInput.m](BuildMeshToComputaionalModel/sWT_1_Male_ComputationalModelInput.m)

---

### 3. Compupational model

Code in: ComputationalModel/

Code for applying preassure to the mesh to deploy it.

Based on FEM and Lagrangian models, this code deploys the quadratic mesh previoulsy built.

Run [Pressure_Expansion_TQQL_Model.m](BuildMeshToComputaionalModel/Pressure_Expansion_TQQL_Model.m) and use as input any ComputationalModelInputs as [sWT_1_Male_ComputationalModelInput.m](BuildMeshToComputaionalModel/sWT_1_Male_ComputationalModelInput.m).

The results will be saved as .mat and .vtk in two different folders.

---

### 4. Extract the adult eyes

Code in: ExtractAdultEyes/

Code for extracting the apical surface of adult eyes.

The segmentation must be load in .tif and the result will be saved in .mat.

Run [ExtractAdultEye.m](BuildMeshToComputaionalModel/ExtractAdultEye.m) to extract the apical surface of any segmented adult eye.

---

### 5. Processing and Analyzing the meshes

Code in: ProcessMeshes/

Code for process the simulations and adult (downsample, smooth and isometric explicit remesh), obtain their Gaussian curvature profiles and their Gaussian metric.

---
