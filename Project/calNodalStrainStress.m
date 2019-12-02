function [strain, stress] = calNodalStrainStress(glU,meshStruct,layer_position)
%nodalGradient = calNodalGradient(glField,meshStruct)   Calculate nodal
%stress and strain
%
% Input
% glU    : column or row vector, the global nodal displacement;
% meshStruct : structure containing FE mesh info, see LinElast.m for
% detailed spec.
%
% Output
% strain: nx3 matrix, first column xx components, second column yy
% components, third column xy components
% stress: nx3 matrix, first column xx components, second column yy
% components, third column xy components

% Algorithm
% This code first calculates stress and strain at gauss points, and
% then extrapolates to the nodal points. That way the output gradients
% are one order more accurate than calculating them directly at nodal
% points. This technique in FE is called "stress recovery", and can be
% applied in both scalar field (2DBVP) and vector field (Linear elastic)
% problems.
% The idea of the extrapolation is global least square. In the end, we form
% M*u=R, where M is the global mass matrix M = integral(N'*N) and
% R = integral(N'*filed_gp). u is the nodal gradient values. See lecture
% notes for more info.

% last edit: Nov 18 2015 Y. Xu

if isrow(glU)
    glU = glU'; % to column vector
end

% unpack things we need
numNodes=meshStruct.numNodes;
nnpe=meshStruct.nnpe;
numEls=meshStruct.numEls;
nCoords=meshStruct.nCoords;
elCon=meshStruct.elCon;
gatherMat=meshStruct.gatherMat;
D=meshStruct.Material.D;
thickness = meshStruct.thickness;
thickness_layer = meshStruct.thickness_layer;

% % Initialization
% R = zeros(numNodes,6);  % right hand side (3 stress components + 3 strain components)
% M = spalloc(numNodes,numNodes,10*numNodes);  % allocate sparse matrix, last input
%                                              % is the estimation of
%                                              % nonzero numbers of M
[qp,w] = Gauss(nnpe);
for e = 1:numEls
    %     % Elemental initialization
    %     Me = zeros(nnpe*3,1);
    %     Re = zeros(nnpe*3,6);
    % get nodal coordinates for this element
    xvec=nCoords(elCon(e,:),1); % these are column vectors of the
    yvec=nCoords(elCon(e,:),2); % physical coordinates of the
    % nodes of this element
    
    x_elemlength =  abs(xvec(2) - xvec(1))/2;
    y_elemlength =  abs(yvec(3) - yvec(2))/2;
    x_elem_center = mean(xvec);
    y_elem_center = mean(yvec);
    glb_nodes = gatherMat(e,:); % local displacement
    
    u_z = d(glb_nodes*3-2);
    theta_x = d(glb_nodes*3-1);
    theta_y = d(glb_nodes*3);
    lcU = [];
    for jj = 1:1:nnpe
        lcU = [lcU; d(glb_nodes(jj)*3-2);d(glb_nodes(jj)*3-1);d(glb_nodes(jj)*3)];
    end
    for iqp = 1:length(w)
        NofXiEta=Nmatrix(qp(iqp,:));  % 1x12 matrix,  row vector of shape functions
        % at this QP in parent domain
        dNdXiEta=dNmatrix(qp(iqp,:)); % 2x12 matrix, of first derivs of shape func
        % at this QP in parent domain
        ddNdXiEta=ddNmatrix(qp(iqp,:)); % 3x12 matrix, of second derivs of shape func
        % at this QP in parent domain
        
        JofXiEta= [x_elemlength, 0; 0, y_elemlength]; % jacobian at this QP (2x2 matrix)
        
        XY= [x_elemlength*qp(iqp,1)+x_elem_center, y_elemlength*qp(iqp,2)+y_elem_center];   % physical coordinate of this QP (1x2)
        
        % now create the N and B matrices
        
        N = [NofXiEta; (1/x_elemlength)*dNdXiEta(1,:); (1/y_elemlength)*dNdXiEta(2,:)];
        
        B =  [(1/x_elemlength)^2*ddNdXiEta(1,:); (1/y_elemlength)^2*ddNdXiEta(2,:); (2/x_elemlength)*(1/y_elemlength)*ddNdXiEta(3,:)];
        % elemental strain and stress
        lcStrain = -thickness_layer(layer_position)*B*lcU;
        lcStress = (D*12/(thickness)^3)* lcStrain;
    end
    
    % Assembly: we can do the same thing as in Assembly.m, but here we
    % use a simple way--a little inefficient for sparse matrix. This would
    % still be better than (when solving the large linear system) using
    % full matrix.
    
    glInd = elCon(e,:); % global index
    M(glInd,glInd) = M(glInd,glInd) + Me;
    R(glInd,:) = R(glInd,:) + Re;
end

strain = M\R(:,1:3);
stress = M\R(:,4:6);



