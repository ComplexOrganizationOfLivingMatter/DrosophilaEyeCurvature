%% MAT-fem_Shells
% 3 Nodes Triangular Thick/Thin Shell Element TQQL 

% Clear memory and variables
  clc
  clear all
    
% The variables are read as a MAT-fem subroutine
% young = Young Modulus
% poiss = Poission Ratio
% thick = Thickness 
% denss = Density
% coordinates = [ x , y , z ] coordinate matrix nnode x ndime (3) 
% elements    = [ inode , jnode , knode ] element connectivity’ matrix.
%               Matrix size: nelem x nnode; nnode = 3
% fixnodes    = [ node number , dof , fixed value ] matrix with 
%               Dirichlet restrictions, were dof=1,2,3 for displacements,
%               dof=4 for rotation in x' and dof=5 for rotation in y'
% pointload   = [ node number , dof , load value ] matrix with
%               nodal loads, were dof=1,2,3 for loads,
%               dof=4 for x' moment and dof=5 for y' moment
% uniload     = [ normal surface load ] sparse matrix size: nelem x 1

  file_name = input('Enter the file name: ','s');

  eval(file_name);       % Read input file

%% First step (initial conditions)

% Find basic dimensions

  npnod = size(coordinates,1);         % Number of nodes
  nelem = size(elements,1);            % Number of elements
  nnode = size(elements,2);            % Number of nodes per element
  dofpn = 5;                           % Number of DOF per node
  dofpe = nnode*dofpn;                 % Number of DOF per element
  nndof = npnod*dofpn;                 % Number of total DOF 

  elements = sortrows(elements);
  
  N_iterations = 1000;

% Dimension the global matrices

  StifMat(:,:,1)  = zeros( nndof , nndof );  % Create the global stiffness matrix
  force(:,1)    = zeros( nndof , 1 );      % Create the global force vector
  force1(:,1)   = zeros( nndof , 1 );      % Create the global force vector
  reaction(:,1) = zeros( nndof , 1 );      % Create the global reaction vector
  u(:,1)        = zeros( nndof , 1 );      % Nodal variables
  
% Material properties (Constant over the domain)

  aux1 = thick*young/(1-poiss^2);
  aux2 = poiss*aux1;
  aux3 = thick*young/2/(1+poiss);
  aux4 = (5/6)*thick*young/2/(1+poiss);
   
  D_matm = [aux1,aux2,   0;
            aux2,aux1,   0;
               0,   0,aux3];
  
  D_mats = [aux4,   0;
               0,aux4];
  
  D_matb = D_matm*(thick^2/12);
  
% Gauss point coordinates

  xg(1) = 1.0/6.0;
  xg(2) = 2.0/3.0;
  xg(3) = 1.0/6.0;

  yg(1) = 1.0/6.0;
  yg(2) = 1.0/6.0;
  yg(3) = 2.0/3.0;

  wg(1) = 1.0/3.0;
  wg(2) = 1.0/3.0;
  wg(3) = 1.0/3.0;

%% Part 2: Global System Assembly

% Element cycle

  for ielem = 1 : nelem

      lnods(1:nnode) = elements(ielem,1:nnode);
      
      % Part 2.1:
  
      cxyz(1:6,1:3) = coordinates(lnods(1:6),1:3);      % Element coordinates
      ctrv = [ cxyz(1,1:3) ;
               cxyz(3,1:3) ;
               cxyz(5,1:3) ];
      Te = Rotation_system_v1_1(ctrv);

      la(ielem,1:3) = Te(1,1:3);
      la(ielem,4:6) = Te(2,1:3);
      la(ielem,7:9) = Te(3,1:3);

      ctxy = cxyz*transpose(Te);    % Rotate coordinates to element mid plane
      ctxy(2,1:2) = (ctxy(1,1:2) + ctxy(3,1:2))/2;   % Assume the middle node
      ctxy(4,1:2) = (ctxy(3,1:2) + ctxy(5,1:2))/2;   % in the same plane
      ctxy(6,1:2) = (ctxy(5,1:2) + ctxy(1,1:2))/2;

      x = ctxy(1:6,1); % Elem. X coordinates
      y = ctxy(1:6,2); % Elem. Y coordinates
     
      K_elem = zeros( dofpe , dofpe );

      areae = 0.0;
      fx = zeros(1,6);
      fy = zeros(1,6);
      fz = zeros(1,6);
      
      % Part 2.2 and 2.3:
    
      for igaus = 1 : 3
          
          [bmat_b,bmat_s,bmat_m,N,area] = B_mat_Shell_TQQL_v1_1(x,y,xg(igaus),yg(igaus),Te); 

          K_b = transpose(bmat_b)*D_matb*bmat_b*area*wg(igaus);
          K_m = transpose(bmat_m)*D_matm*bmat_m*area*wg(igaus);
          K_s = transpose(bmat_s)*D_mats*bmat_s*area*wg(igaus);
      
          K_elem = K_elem + K_b + K_m + K_s;
      
          areae = areae + area*wg(igaus);
          fx = fx + uniload(ielem)*area*Te(3,1)*wg(igaus)*N;
          fy = fy + uniload(ielem)*area*Te(3,2)*wg(igaus)*N;
          fz = fz + (uniload(ielem)*Te(3,3) - denss*thick)*area*wg(igaus)*N;

      end
    
      ElemFor = [fx(1),fy(1),fz(1),0,0,fx(2),fy(2),fz(2),0,0,...
                 fx(3),fy(3),fz(3),0,0,fx(4),fy(4),fz(4),0,0,...
                 fx(5),fy(5),fz(5),0,0,fx(6),fy(6),fz(6),0,0];

% Find the equation number list for the i-th element

      for i = 1 : nnode

          ii = (i-1)*dofpn;

          for j = 1:dofpn

              eqnum(ii+j) = (lnods(i)-1)*dofpn + j;   % Build the eq. number list

          end
      end

% Part 2.4:

% Assemble the force vector and the stiffness matrix

      for i = 1 : dofpe

          ipos = eqnum(i);
          force(ipos,1) = force(ipos,1) + ElemFor(i);

          for j = 1 : dofpe

              jpos = eqnum(j);
              StifMat (ipos,jpos,1) = StifMat (ipos,jpos,1) + K_elem(i,j);

          end
      end

  end  % End element cycle

%% Part 3: Viscous Regularization

% Thanks to Shell_TQQL_v1_1_iterations.m we found out that our system has
% an inestability, this means that in our system there is a movement as
% a rigid solid and our matrix system (StifMat * u = force) is not
% solved correctly. As a result, our triangles start to rotate so many
% times.

% In order to solve this inestability we will introduce in our matrix
% system a viscous regularization term, which physically means that we
% will introduce our triangles into a viscous medium. In this way, they
% will experience a force, which is proportional to the speed of the
% triangles. This force will not allow them to rotate so much and will
% solve the inestability.
% With all of this, the matrix system will have the following form:
%
%                (c+StifMat)*u_t+1 = force_t+1 + c*u_t (1)
% 
% where c:= viscous regularization parameter and will have the following
% form:
%
%                    c= alpha*OnesDisplacements
% 
% where alpha:= parameter that need to be calibrate (as a starting point
% it will be the elements of StifMat that refer to the displacement,
% which are 1,2,3,6,7,8....,*10^(-3)). It's obvious that alpha needs to
% be small because it's a perturbation in the problem and we don't want
% to really change the problem, we just want to make it coherent.
%
% where OnesDisplacement:= it is a diagonal matrix nndof x nndof that
% will have 1 is the displacements positions and 0 in other case. This
% means that elements 1,2,3,6,7,8,11,12,13... in the diagonal will be 1
% and the rest that is out of this progression will be 0.

% We will start with OnesDisplacements.

  Displacements_diagonal = zeros(1,nndof);

  Displacementes_diagonal(1:3) = 1;

  for i = 0:(npnod-1)

      Displacements_diagonal(5*i+1:5*i+3) = 1;

  end

% So the matrix will be:

  OnesDisplacements = diag(Displacements_diagonal); 

% And alpha will be the medium value of StifMat*10^(-3).

  alpha = (max(max(full(StifMat(:,:,1))))+min(min(full(StifMat(:,:,1)))))/2*10^(-6);

% So the viscous regularization parameter is:

  c = alpha*OnesDisplacements;

%% Part 4: Application of Dirichlet Boundary Conditions and Load Vector Assembly

% Part 4.1:

% Add point load conditions to the force vector

  for i = 1 : size(pointload,1)

      ieqn = (pointload(i,1)-1)*dofpn + pointload(i,2);   % Find eq. number
      force(ieqn,1) = force(ieqn,1) + pointload(i,3);     % and add the force

  end
 
% Part 4.2:

% Apply the Dirichlet conditions and adjust the right hand side

  j = 0;

  for i = 1 : size(fixnodes,1)

      ieqn = (fixnodes(i,1)-1)*dofpn + fixnodes(i,2);  % Find equation number
      u(ieqn,1) = fixnodes(i,3);                  % and store the solution in u
      j = j + 1;
      fix(j) = ieqn;                        % and mark the eq. as a fix value

  end
  
  force1(:,1) = force(:,1) - (c + StifMat(:,:,1)) * u(:,1);      % Adjust the rhs with the known values

%% Part 5: Linear System Solution

% Now we need to adapt the matrix system to the new scheme. In the previous
% code we did:
%
% Compute the solution by solving StifMat * u = force for the remaining
% unknown values of u
%
%  u(FreeNodes,w) = StifMat(FreeNodes,FreeNodes,w) \ force1(FreeNodes,w);
%
% Now our system will be the one mentioned above (1) so, using the
% nomenclature the matrix system will be:
%
%        (c+StifMat)*u(:,w) = force1(:,w) + c*u(:,w-1)  (2)
%


% We will now solve this matrix system.

  FreeNodes = setdiff( 1:nndof , fix );           % Find the free node list
                                                  % and solve for it

  NewStifMat = c + StifMat(:,:,1);

% So, the new nodal variables are:

  u(FreeNodes,1) = NewStifMat(FreeNodes,FreeNodes)\force1(FreeNodes,1);

% Compute the stresses

  [Strnod,GaussPoints_Stresses] = Stress_Shell_TQQL_v1_1(D_matb,D_matm,D_mats,xg,yg,u(:,1),coordinates,elements);

%% The calculation of the membrane enery won't be taken into account for the
% paper

  aux_parameter_1 = 1/(1+poiss)^2;
  aux_parameter_2 = 1/(1-poiss^2);

% For convenience:

  for i = 1:size(GaussPoints_Stresses,2)

      GaussPoints_Stresses{i}(:,4:6) = GaussPoints_Stresses{i}(:,4:6)/aux1;

  end
%% Area Calculation

% For the area we will use the Heron's formula that doesn't need the height
% but instead it uses the semiperimeter. We will create elem_vert (elements
% vertices) to make easier the area calculus:

  elem_vert = [elements(:,1),elements(:,3),elements(:,5)];
  Trian_Area = zeros(N_iterations,nelem);

  for i = 1:nelem

      for pg = 1:3
          
          Membrane_Energy_Gauss(pg) = (1/2)*(aux_parameter_1*(GaussPoints_Stresses{i}(1,pg+3)+GaussPoints_Stresses{i}(2,pg+3))^2 - 2*(1-poiss)*aux_parameter_2*((GaussPoints_Stresses{i}(1,pg+3)-poiss*GaussPoints_Stresses{i}(2,pg+3))*(GaussPoints_Stresses{i}(2,pg+3)-poiss*GaussPoints_Stresses{i}(1,pg+3))-GaussPoints_Stresses{i}(3,pg+3)^2))*wg(pg);
      
      end

      Membrane_Energy_ite(1,i) = sum(Membrane_Energy_Gauss);
      Membrane_Energy_tot = sum(Membrane_Energy_ite);           
      Trian_Area(1,i) = Triangle_Area_NonLinearTriangulation(elements(i,:),coordinates);
    
  end
  
%% Part 1: Normal Directions Update and Load Update Procedure

% I considerer the begin of a iteration the calculation of new normals even though it appears at the end.

% We need to update the coordinates. Every node has 5 degrees of freedom,
% that's why u is a (nnode*5)x1 matrix, and the first 3 are the
% displacements in x,y,z respectively; so we will create a matrix with
% these componentes, for example, the first row will be 1,2,3; the second
% row will be 6,7,8 and so on.

% The first thing we need to do is access to the elements of u, which is a
% sparse matrix.

  u_aux = full(u(:,1));

% Now, u_aux contains the elements of u. For example u(3) = (3,1) 0.5, so
% we have that u_aux(3) = 0.5;

% Now we will transform u_aux in a matrix npnod x 3, where every row has
% the componentes 1,2,3; 6,7,8;11,12,13... and so on, of u.

% The first 3 will be outside of the bucle (it's easier).

  displacements(1,:,1) = u_aux(1:3)';

  for k=0:(npnod-1)

      displacements(k+1,:,1) = u_aux(dofpn*k+1:dofpn*k+3)';

  end
  
%% Save Results and Graphic representation:

% We need the coordinates of the vetices, that are:

  path2save = pwd;
  filename = strcat(file_name,'_Pressure_Area_Original_',num2str(1));
  mkdir(strcat(pwd,'/','Results/Pressure/Paraview_Results/Area_Original_t_0125_p_1e5/'))
  mkdir(strcat(pwd,'/','Results/Pressure/Paraview_Results/Area_Evolution_t_0125_p_1e5/'))
  mkdir(strcat(pwd,'/Results/Pressure/MATLAB_Results/'))
  
%   vtkfilewriter(strcat(pwd,'/','Results/Pressure/Paraview_Results/Area_Original_t_0125_p_1e5/',filename),'ASCII',coordinates,elem_vert,Trian_Area(1,:),Membrane_Energy_tot);
%   filename_Area = strcat(file_name,'_Pressure_Area_Evolution_',num2str(1));
%   vtkfilewriter(strcat(pwd,'/Results/Pressure/Paraview_Results/Area_Evolution_t_0125_p_1e5/',filename_Area),'ASCII',coordinates,elem_vert,Trian_Area(1,:),Membrane_Energy_tot);
  save(strcat(pwd,'/Results/Pressure/MATLAB_Results/',file_name,'_Pressure_',num2str(1),'_Variables.mat'),'coordinates','elements','-v7.3');
  
%% Still part 1:
  
  % Now we summ coordinates and the displacements.

  coordinates = coordinates + displacements(:,:,1);

  % Now, we obtain the normals:

  elementsVertices = [elements(:,1),elements(:,3),elements(:,5)];
  vertices = unique(elementsVertices);
  vertices = setdiff(vertices,unique(fixnodes(:,1)));
  coordinatesVertices = coordinates(vertices,:);
  normals = findNormals(coordinatesVertices);

  % Now we adjust the load to be like a preassure. The original load will
  % be:

  originalLoad = pointload(1,3);

  % So, pointload for the next iteration will be:

  pointload = FixPressurePoinload(pointload,normals,originalLoad);


%% Iterations

% Now we will solve for \Delta u because the equations are easier to solve
% numerically), so now u will be Delta_u (the equations need to be changed)
% and uu will be u. In order to updarte the coordinates, we will use that
% we know the following:
%
%                     uu(:,w) = uu(:,w-1) + u(:,w)
%
% So. now, uu = u and u will be \Delta u: 

  uu = u; 

  for w = 2:N_iterations

      % Find basic dimensions

      npnod = size(coordinates,1);         % Number of nodes
      nelem = size(elements,1);            % Number of elements
      nnode = size(elements,2);            % Number of nodes per element
      dofpn = 5;                           % Number of DOF per node
      dofpe = nnode*dofpn;                 % Number of DOF per element
      nndof = npnod*dofpn;                 % Number of total DOF

      StifMat  = zeros( nndof , nndof );  % Create the global stiffness matrix
      force    = zeros( nndof , 1 );      % Create the global force vector
      force1   = zeros( nndof , 1 );      % Create the global force vector
      reaction = zeros( nndof , 1 );      % Create the global reaction vector
      u(:,w)        = zeros( nndof , 1 );      % Nodal variables

      elements = sortrows(elements);
      
  % Element cycle

      for ielem = 1 : nelem
          
          lnods(1:nnode) = elements(ielem,1:nnode);
  
          cxyz(1:6,1:3) = coordinates(lnods(1:6),1:3);      % Element coordinates
          ctrv = [ cxyz(1,1:3) ;
                   cxyz(3,1:3) ;
                  cxyz(5,1:3) ];
          Te = Rotation_system_v1_1(ctrv);

          la(ielem,1:3) = Te(1,1:3);
          la(ielem,4:6) = Te(2,1:3);
          la(ielem,7:9) = Te(3,1:3);

          ctxy = cxyz*transpose(Te);    % Rotate coordinates to element mid plane
          ctxy(2,1:2) = (ctxy(1,1:2) + ctxy(3,1:2))/2;   % Assume the middle node
          ctxy(4,1:2) = (ctxy(3,1:2) + ctxy(5,1:2))/2;   % in the same plane
          ctxy(6,1:2) = (ctxy(5,1:2) + ctxy(1,1:2))/2;

          x = ctxy(1:6,1); % Elem. X coordinates
          y = ctxy(1:6,2); % Elem. Y coordinates
     
          K_elem = zeros( dofpe , dofpe );

          areae = 0.0;
          fx = zeros(1,6);
          fy = zeros(1,6);
          fz = zeros(1,6);
    
          for igaus = 1 : 3

              [bmat_b,bmat_s,bmat_m,N,area] = B_mat_Shell_TQQL_v1_1(x,y,xg(igaus),yg(igaus),Te); 

              K_b = transpose(bmat_b)*D_matb*bmat_b*area*wg(igaus);
              K_m = transpose(bmat_m)*D_matm*bmat_m*area*wg(igaus);
              K_s = transpose(bmat_s)*D_mats*bmat_s*area*wg(igaus);
      
              K_elem = K_elem + K_b + K_m + K_s;
      
              areae = areae + area*wg(igaus);
              fx = fx + uniload(ielem)*area*Te(3,1)*wg(igaus)*N;
              fy = fy + uniload(ielem)*area*Te(3,2)*wg(igaus)*N;
              fz = fz + (uniload(ielem)*Te(3,3) - denss*thick)*area*wg(igaus)*N;

          end
    
          ElemFor = [fx(1),fy(1),fz(1),0,0,fx(2),fy(2),fz(2),0,0,...
                     fx(3),fy(3),fz(3),0,0,fx(4),fy(4),fz(4),0,0,...
                     fx(5),fy(5),fz(5),0,0,fx(6),fy(6),fz(6),0,0];

% Find the equation number list for the i-th element

          for i = 1 : nnode

              ii = (i-1)*dofpn;

              for j =1:dofpn

                  eqnum(ii+j) = (lnods(i)-1)*dofpn + j;   % Build the eq. number list
                  
              end
          end
 
% Assemble the force vector and the stiffness matrix

          for i = 1 : dofpe

              ipos = eqnum(i);
              force(ipos) = force(ipos) + ElemFor(i);

              for j = 1 : dofpe
                  
                  jpos = eqnum(j);
                  StifMat (ipos,jpos) = StifMat (ipos,jpos) + K_elem(i,j);
                  
              end
          end

      end  % End element cycle

% Add point load conditions to the force vector

      for i = 1 : size(pointload,1)

          ieqn = (pointload(i,1)-1)*dofpn + pointload(i,2);   % Find eq. number
          force(ieqn) = force(ieqn) + pointload(i,3);     % and add the force

      end

% Apply the Dirichlet conditions and adjust the right hand side

      j = 0;

      for i = 1 : size(fixnodes,1)

          ieqn = (fixnodes(i,1)-1)*dofpn + fixnodes(i,2);  % Find equation number
          u(ieqn,w) = fixnodes(i,3);                  % and store the solution in u
          j = j + 1;
          fix(j) = ieqn;                        % and mark the eq. as a fix value

      end
  
      force1(:) = force(:) - (c + StifMat) * u(:,w) - StifMat*uu(:,w-1);      % Adjust the rhs with the known values

% Now we need to adapt the matrix system to the new scheme. In the previous
% code we did:
%
% Compute the solution by solving StifMat * u = force for the remaining
% unknown values of u
%
%  u(FreeNodes,w) = StifMat(FreeNodes,FreeNodes,w) \ force1(FreeNodes,w);
%
% Now our system will be the one mentioned above (1) so, using the
% nomenclature the matrix system will be:
% coordinates(elements(1,1),:)
%        (c+StifMat)*u(:,w) = force1(:,w) + c*u(:,w-1)  (2)
%
% We will now solve this matrix system.

      FreeNodes = setdiff( 1:nndof , fix );           % Find the free node list
                                                  % and solve for it

      NewStifMat = c + StifMat;

% So, the new nodal variables are:

      u(FreeNodes,w) = NewStifMat(FreeNodes,FreeNodes)\(force1(FreeNodes));

% So the change in nodal variables will be:

% Delta_u = u(:,w);
  
      uu(:,w) = uu(:,w-1) + u(:,w);
  
      Delta_u = uu(:,w);

% We need the full nodal variables so:

      full_nodal_variables = sum(uu(:,1:w),2); 

% Compute the stresses

      [Strnod,GaussPoints_Stresses] = Stress_Shell_TQQL_v1_1(D_matb,D_matm,D_mats,xg,yg,Delta_u,coordinates,elements);

% For convenience:

      for i = 1:size(GaussPoints_Stresses,2)

          GaussPoints_Stresses{i}(:,4:6) = GaussPoints_Stresses{i}(:,4:6)/aux1;

      end

% For the area we will use the Heron's formula that doesn't need the height
% but instead it uses the semiperimeter. We will create elem_vert (elements
% vertices) to make easier the area calculus:

      elem_vert = [elements(:,1),elements(:,3),elements(:,5)];
    
      for i = 1:nelem

          Trian_Area(w,i) = Triangle_Area_NonLinearTriangulation(elements(i,:),coordinates);
    
       end

 
% We need to update the coordinates. Every node has 5 degrees of freedom,
% that's why u is a (nnode*5)x1 matrix, and the first 3 are the
% displacements in x,y,z respectively; so we will create a matrix with
% these componentes, for example, the first row will be 1,2,3; the second
% row will be 6,7,8 and so on.

% The first thing we need to do is access to the elements of u, which is a
% sparse matrix.

      u_aux = full(Delta_u); 

% Now, u_aux contains the elements of u. For example u(3) = (3,1) 0.5, so
% we have that u_aux(3) = 0.5;

% Now we will transform u_aux in a matrix npnod x 3, where every row has
% the componentes 1,2,3; 6,7,8;11,12,13... and so on, of u.

% The first 3 will be outside of the bucle (it's easier).

      displacements(1,:,w) = u_aux(1:3)';
    
      for k=0:(npnod-1)

          displacements(k+1,:,w) = u_aux(dofpn*k+1:dofpn*k+3)';
    
      end

% Graphic representation

  filename = strcat(file_name,'_Pressure_Area_Original_',num2str(w));

%   vtkfilewriter(strcat(pwd,'/','Results/Pressure/Paraview_Results/Area_Original_t_0125_p_1e5/',filename),'ASCII',coordinates,elem_vert,Trian_Area(1,:),Membrane_Energy_tot);
% 
%   filename_Area = strcat(file_name,'_Pressure_Area_Evolution_',num2str(w));
% 
%   vtkfilewriter(strcat(pwd,'/Results/Pressure/Paraview_Results/Area_Evolution_t_0125_p_1e5/',filename_Area),'ASCII',coordinates,elem_vert,Trian_Area(w,:),Membrane_Energy_tot);

  save(strcat(pwd,'/Results/Pressure/MATLAB_Results/',file_name,'_Pressure_',num2str(w),'_Variables.mat'),'coordinates','elements','-v7.3');

  % ToGiD_Shell_TQQL_v1_1(file_name,u(:,1),reaction(:,1),Strnod,1,initial_coordinates,elements,displacements,GaussPoints_Stresses); 
  
  % Now we summ coordinates and the displacements.

  coordinates = coordinates + displacements(:,:,w);

  % Now, we obtain the normals:

  coordinatesVertices = coordinates(vertices,:);
  normals = findNormals(coordinatesVertices);

  % The pointload matrix is already well-built so, we just need to change
  % the third column with the new normals:

  normalsNew = originalLoad*reshape(normals.',1,[])';

  pointload(:,3) = normalsNew;

  end

  save(strcat(file_name,'_Variables.mat'),'coordinates','elements','Trian_Area','Membrane_Energy_tot','-v7.3')
