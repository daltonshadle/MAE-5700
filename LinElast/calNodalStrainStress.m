function [strain, stress] = calNodalStrainStress(glU,meshStruct)
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

% Initialization
R = zeros(numNodes,6);  % right hand side (3 stress components + 3 strain components)
M = spalloc(numNodes,numNodes,10*numNodes);  % allocate sparse matrix, last input
                                             % is the estimation of
                                             % nonzero numbers of M
[qp,w] = Gauss(nnpe);
for e = 1:numEls
    % Elemental initialization
    Me = zeros(nnpe);        
    Re = zeros(nnpe,6);
    
    xynode = nCoords(elCon(e,:),:); % nodal coordinate
    lcU = glU(gatherMat(e,:)); % local displacement
    for iqp = 1:length(w)
        NofXiEta=Nmatrix(qp(iqp,:),nnpe); % row vector of shape functions
                                          % at this QP in parent domain
        dNdXiEta=dNmatrix(qp(iqp,:),nnpe);% 2xnnpe array of shape func derivs
                                      % at this QP in parent domain
        JofXiEta=dNdXiEta*xynode;    % jacobian at this QP (2x2 matrix)
        dNdXY=JofXiEta\dNdXiEta;     % 2xnnpe array of dNdX at ths QP
     
        % Here we no longer need the N matrix as in TwoDElem.m
        B=[];
        for np=1:nnpe         
            B=[B,[dNdXY(1,np), 0; 0, dNdXY(2,np); dNdXY(2,np), dNdXY(1,np)]];
        end
        % elemental strain and stress
        lcStrain = B * lcU;
        lcStress = D * lcStrain;
        
        Me = Me + NofXiEta'*NofXiEta * w(iqp)*det(JofXiEta);
        Re = Re + NofXiEta'* [lcStrain',lcStress'] * w(iqp)*det(JofXiEta);

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



