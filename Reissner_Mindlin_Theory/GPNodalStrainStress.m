function [strain, stress, moment] = GPNodalStrainStress(glU,meshStruct,layer_position)
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
Db=meshStruct.Material.Db;
Ds=meshStruct.Material.Ds;
sc=meshStruct.Material.sc;
thickness=meshStruct.thickness;

% Initialization
R = sparse(numNodes,13);  % right hand side (5 stress components + 5 strain components + 3 moments)
M = sparse(numNodes,numNodes);

[qp,w] = Gauss(nnpe);
for e = 1:numEls
    % Elemental initialization
    Me = zeros(nnpe);
    Re = zeros(nnpe,13);
    
    % get nodal coordinates for this element
    xvec=nCoords(elCon(e,:),1); % these are column vectors of the
    yvec=nCoords(elCon(e,:),2); % physical coordinates of the
    % nodes of this element
    
    
    
    lcU = glU(gatherMat(e,:)); % local displacement
    for iqp = 1:length(w)
        % for each quadrature point ...
        NofXiEta=Nmatrix(qp(iqp,:)); % row vector of shape functions
        % at this QP in parent domain (1x4)
        dNdXiEta=dNmatrix(qp(iqp,:));% 2xnnpe array of shape func derivs
        % at this QP in parent domain (2x4)
        JofXiEta=dNdXiEta*[xvec,yvec];    % jacobian at this QP (2x2 matrix)
        dNdXY=inv(JofXiEta)*dNdXiEta;     % 2xnnpe array of dNdX at ths QP
        
        XY=NofXiEta*[xvec,yvec]; % physical coordinate of this QP (1x2)
        
        % now create the N and Bb (bending) and Bs (shear) matrices
        N=[]; Bb=[]; Bs=[];
        for np=1:nnpe
            N=[N,[NofXiEta(np), 0, 0; 0, NofXiEta(np), 0; 0, 0, NofXiEta(np)]];
            Bb=[Bb,[0, dNdXY(1,np), 0; 0, 0, dNdXY(2,np); 0, dNdXY(2,np), dNdXY(1,np)]];
            Bs=[Bs,[dNdXY(1,np), -NofXiEta(np), 0; dNdXY(2,np), 0, -NofXiEta(np)]];
        end
        
        % elemental strain and stress
        lcStrain_b = - layer_position * Bb * lcU;
        lcStress_b = Db * lcStrain_b;
        lcMoment = - Db * thickness^3 * Bb * lcU / 12;
        
        lcStrain_s = Bs * lcU;
        lcStress_s = Ds * lcStrain_s;
        
        lcStrain = [lcStrain_b; lcStrain_s];
        lcStress = [lcStress_b; lcStress_s];
        
        Me = Me + NofXiEta'*NofXiEta * w(iqp)*det(JofXiEta);
        Re = Re + NofXiEta'* [lcStrain',lcStress',lcMoment'] * w(iqp)*det(JofXiEta);
        
    end
    
    % Assembly: we can do the same thing as in Assembly.m, but here we
    % use a simple way--a little inefficient for sparse matrix. This would
    % still be better than (when solving the large linear system) using
    % full matrix.
    
    glInd = elCon(e,:); % global index
    M(glInd,glInd) = M(glInd,glInd) + Me;
    R(glInd,:) = R(glInd,:) + Re;
end

strain = M\R(:,1:5);
stress = M\R(:,6:10);
moment = M\R(:,11:13);
end


