function [strain, stress] = AtNodalStrainStress(glU,meshStruct,layer_position)
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
thickness=meshStruct.thickness

% Initialization
strain = zeros(numNodes,3);
stress = zeros(numNodes,3);

for e = 1:numEls
    % Elemental initialization

    % get nodal coordinates for this element
    xvec=nCoords(elCon(e,:),1); % these are column vectors of the
    yvec=nCoords(elCon(e,:),2); % physical coordinates of the
                                    % nodes of this element

    x_elemlength =  abs(xvec(2) - xvec(1))/2;                           
    y_elemlength =  abs(yvec(3) - yvec(2))/2;         
    x_elem_center = mean(xvec);
    y_elem_center = mean(yvec);
    
    
    lcU = glU(gatherMat(e,:)); % local displacement
    
    pos = [(xvec - x_elem_center)/x_elemlength, (yvec - y_elem_center)/y_elemlength];
    
    lcStrain = zeros(3,4);
    lcStress = zeros(3,4);
    
    for i = 1:size(pos,1)
        % Calcluate ddN matrix
        ddNdXiEta=ddNmatrix(pos(i,:)); % 3x12 matrix, of second derivs of shape func
                                        % at this QP in parent domain

        % now create the B matrix
        B =  [(1/x_elemlength)^2*ddNdXiEta(1,:); (1/y_elemlength)^2*ddNdXiEta(2,:); ...
            (2/x_elemlength)*(1/y_elemlength)*ddNdXiEta(3,:)];

        % elemental strain and stress 3x4 matrices
        lcStrain(:,i) = - layer_position * B * lcU;
        lcStress(:,i) = 12 * D * lcStrain(:,i) / thickness^3;

    end
    lcStrain = lcStrain';
    lcStress = lcStress';
    glInd = elCon(e,:); % global index

    strain(glInd,:) = strain(glInd,:) + lcStrain;
    stress(glInd,:) = stress(glInd,:) + lcStress;
end
end



