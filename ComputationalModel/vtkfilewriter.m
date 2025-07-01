function [vtkfile] = vtkfilewriter(filename,format,coordinates,elements,Trian_Area,Membrane_Energy)

%% Inputs

% filename = string. It's the name of the .vtk file.

% format = the format to save the .vtk file. It can be 'ASCII' or 'BINARY'.

% coordinates = matrix of points that conforms the mesh. (Nx3 where N is
%               the number of points).

% elements = connectivity matrix. Its a Mx3 matrix (M = number of
%            triangles) where row i has the indices of the rows of coordinates
%            that are the vertices of the triangle i.

% Trian_Area = Area of the triangles. It's a vector whose i component is
%              the area of the triangle i in the triangulation.

% Membrane_Energy = Membrane energy of the triangles. It's a vector whose i
%                   component is the membrane energy of the triangle i in
%                   the triangulation.

%% Outputs

% filename.vtk = file in .vtk format for paraview.

%% Code

% We will generate the .vtk file for paraview.

% First, we open the file and give writing permission ('w'):

filename_aux = strcat(filename,'.vtk');

fid = fopen(filename_aux,'w');

% We now write the header and the title:

fprintf(fid,'# vtk DataFile Version 2.0 \n');
filename_aux = strcat(filename,'\n');
fprintf(fid,filename_aux);

% Now the format:

filename_aux = strcat(format,'\n');

fprintf(fid,filename_aux);

% Let's keep on with dataset type, which will be POLYDATA:

fprintf(fid,'DATASET POLYDATA \n');

% We know give the points. To be understood as points we need to identify
% them as POINTS N (number of points) float (the format):

N = size(coordinates,1);

fprintf(fid,'POINTS %d float \n',N);

% And now the matrix of coordinates:

for i = 1:N
    fprintf(fid,'%6.3f %6.3f %6.3f \n',coordinates(i,1),coordinates(i,2),coordinates(i,3));
end

% Now we need to introduce how these points are linked, shaping the
% triangles. As the index starts at 0 we need to substract 1 from elements:

elements = elements-1;

fprintf(fid,'\n');

M = size(elements,1);

fprintf(fid,'POLYGONS %d %d \n',M,4*M);

for i = 1:M
    fprintf(fid,'3 %d %d %d \n',elements(i,1),elements(i,2),elements(i,3));
end

fprintf(fid,'\n');

% For the color of the triangles (area) we need the cell information:

M_1 = size(Trian_Area,2);

fprintf(fid,'CELL_DATA %d \n',M_1);
fprintf(fid,'SCALARS Trian_Area float \n');
fprintf(fid,'LOOKUP_TABLE my_table \n');
for i = 1:M_1
    fprintf(fid,'%0.4f \n',Trian_Area(i));
end

fprintf(fid,'\n');

fprintf(fid,'SCALARS Membrane float \n');
fprintf(fid,'LOOKUP_TABLE my_table_1 \n');
for i = 1:M_1
    fprintf(fid,'%d \n',Membrane_Energy(i));
end

fclose(fid);

