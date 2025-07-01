function Te = Rotation_system_v1_1(cxyz)

%% Rotation_system Evaluates the coordinate system for each element
%
%  Parameters:
%
%    Input, cxyz : Contains the 3D coordinates of the nodes
%                  cxyz(node,axe) 
%   
%    Output, Te the rotation matrix

%
%****  ELEMENT
%

  v12(1) = cxyz(2,1) - cxyz(1,1);
  v12(2) = cxyz(2,2) - cxyz(1,2);
  v12(3) = cxyz(2,3) - cxyz(1,3);

  v13(1) = cxyz(3,1) - cxyz(1,1);
  v13(2) = cxyz(3,2) - cxyz(1,2);
  v13(3) = cxyz(3,3) - cxyz(1,3);

  vze(1) = v12(2)*v13(3) - v12(3)*v13(2);
  vze(2) = v12(3)*v13(1) - v12(1)*v13(3);
  vze(3) = v12(1)*v13(2) - v12(2)*v13(1);
  
  dz = sqrt(vze(1)^2 + vze(2)^2 + vze(3)^2);
  
% Unit vector normal to element surface
  vze(1) = vze(1) / dz;
  vze(2) = vze(2) / dz;
  vze(3) = vze(3) / dz;
  
% XZ plane intesection with element surface
  vxe(1) =  1/sqrt(1+(vze(1)/vze(3))^2); 
  vxe(2) =  0; 
  vxe(3) = -1/sqrt(1+(vze(3)/vze(1))^2); 
  
  dd = vxe(1)*vze(1) + vxe(3)*vze(3);
  if (abs(dd) > 1e-8) 
    vxe(3) = -vxe(3);
  end
  
  if ((vze(3) == 0) && (vze(1) == 0)) 
    vxe(1) =  1; 
    vxe(2) =  0; 
    vxe(3) =  0; 
  end
  
% Vector product
  vye(1) = vze(2)*vxe(3) - vxe(2)*vze(3);
  vye(2) = vze(3)*vxe(1) - vxe(3)*vze(1);
  vye(3) = vze(1)*vxe(2) - vxe(1)*vze(2);
  
  dy = sqrt(vye(1)^2 + vye(2)^2 + vye(3)^2);
  vye(1) = vye(1) / dy;
  vye(2) = vye(2) / dy;
  vye(3) = vye(3) / dy;
  
  if (vye(2) < 0 )
    vye(1) = -vye(1);
    vye(2) = -vye(2);
    vye(3) = -vye(3);
    vxe(1) = -vxe(1);
    vxe(2) = -vxe(2);
    vxe(3) = -vxe(3);
  end
  
  Te = [ vxe(1), vxe(2), vxe(3) ;
         vye(1), vye(2), vye(3) ;
         vze(1), vze(2), vze(3) ];

     
