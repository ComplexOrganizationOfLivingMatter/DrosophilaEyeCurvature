function [Strnod,GaussPoints_Stresses] = Stress_Shell_TQQL_v1_1(D_matb,D_matm,D_mats,xg,yg,u,coordinates,elements)

%% Stress Evaluates the stresses at the Gauss points and smooth the values
%         to the nodes
%
%  Parameters:
%
%    Input, D_matb : Constitutive matrix for bending moment
%           D_matm : Constitutive matrix for membrane forces
%           D_mats : Constitutive matrix for shear force
%           xg     : Local X coordinates of the Gauss points
%           yg     : Local Y coordinates of the Gauss points
%           u      : Nodal displacements
%   
%    Output, Strnod the nodal stress matrix (nnode,nstrs)
%            GaussPoints_Stresses: every cell is the bending moments (first
%                                  3 columns) and membrane momments (last 3
%                                  columns).
  

% Find basic dimensions
  nelem  = size(elements,1);          % Number of elements
  nnode  = size(elements,2);          % Number of nodes por element
  npnod  = size(coordinates,1);       % Number of nodes
  Strnod = zeros(npnod,9);            % Create array for stresses
  dofpn  = 5;                         % Number of DOF per node
  dofpe  = dofpn*nnode;               % Number of DOF per element
  eqnum  = zeros(dofpe);              % Equation number list

% Element cycle

  g_or = [1,3,5]; % vertex position in elements
  for ielem = 1 : nelem
 
    lnods(1:nnode) = elements(ielem,1:nnode);
    
% Find the equation number list for the i-th element
    for i = 1 : nnode
      ii = (i-1)*dofpn;
      for j =1:dofpn
        eqnum(ii+j) = (lnods(i)-1)*dofpn + j;   % Build the eq. number list
      end
    end
    
% Recover the nodal displacements for the i-th element
    u_elem(1:dofpe) = u(eqnum(1:dofpe));
  
    lnods(1:nnode) = elements(ielem,1:nnode);
  
    cxyz(1:6,1:3) = coordinates(lnods(1:6),1:3);      % Element coordinates
    ctrv = [ cxyz(1,1:3) ;
             cxyz(3,1:3) ;
             cxyz(5,1:3) ];
    Te = Rotation_system_v1_1(ctrv);

    ctxy = cxyz*transpose(Te);    % Rotate coordinates to element mid plane
    ctxy(2,1:2) = (ctxy(1,1:2) + ctxy(3,1:2))/2;   % Assume the middle node
    ctxy(4,1:2) = (ctxy(3,1:2) + ctxy(5,1:2))/2;   % in the same plane
    ctxy(6,1:2) = (ctxy(5,1:2) + ctxy(1,1:2))/2;

    x = ctxy(1:6,1);                                  % Elem. X coordinates
    y = ctxy(1:6,2);                                  % Elem. Y coordinates
    
    for igaus = 1 : 3
 
      [bmat_b,bmat_s,bmat_m,N,area] = B_mat_Shell_TQQL_v1_1(x,y,xg(igaus),yg(igaus),Te); 

      Str1(:,igaus)=D_matb*bmat_b*transpose(u_elem);
      Str2(:,igaus)=D_mats*bmat_s*transpose(u_elem);
      Str3(:,igaus)=D_matm*bmat_m*transpose(u_elem);

      % These are the stresses in the gauss points. The code assigns them 
      % one node (without interpolation, that is the theoretically correct)
      % and then we summ over all the triangles that share the same node.
      % So we will store Str1 (bending M) and Str3 (membrane N) for all the
      % gauss points that belong to a triangle.      
      
      Strnod(lnods(g_or(igaus)),1) = Strnod(lnods(g_or(igaus)),1) + Str1(1,igaus);
      Strnod(lnods(g_or(igaus)),2) = Strnod(lnods(g_or(igaus)),2) + Str1(2,igaus);
      Strnod(lnods(g_or(igaus)),3) = Strnod(lnods(g_or(igaus)),3) + Str1(3,igaus);
      Strnod(lnods(g_or(igaus)),4) = Strnod(lnods(g_or(igaus)),4) + Str3(1,igaus);
      Strnod(lnods(g_or(igaus)),5) = Strnod(lnods(g_or(igaus)),5) + Str3(2,igaus);
      Strnod(lnods(g_or(igaus)),6) = Strnod(lnods(g_or(igaus)),6) + Str3(3,igaus);
      Strnod(lnods(g_or(igaus)),7) = Strnod(lnods(g_or(igaus)),7) + Str2(1,igaus);
      Strnod(lnods(g_or(igaus)),8) = Strnod(lnods(g_or(igaus)),8) + Str2(2,igaus);
      Strnod(lnods(g_or(igaus)),9) = Strnod(lnods(g_or(igaus)),9) + 1;
 
    end
    GaussPoints_Stresses{ielem} = [Str1,Str3];

  end
  
  for i = 1 : npnod
    if (Strnod(i,9) ~= 0.0)
      Strnod(i,1:8) = Strnod(i,1:8)/Strnod(i,9);
    end
  end
  
  for ielem = 1 : nelem
    lnods(1:nnode) = elements(ielem,1:nnode);
    Strnod(lnods(2),1:8) = (Strnod(lnods(1),1:8) + Strnod(lnods(3),1:8))/2;
    Strnod(lnods(4),1:8) = (Strnod(lnods(3),1:8) + Strnod(lnods(5),1:8))/2;
    Strnod(lnods(6),1:8) = (Strnod(lnods(5),1:8) + Strnod(lnods(1),1:8))/2;
  end
  
