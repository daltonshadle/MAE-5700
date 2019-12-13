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
D=meshStruct.Material.D;
thickness=meshStruct.thickness;

% Initialization
R = sparse(numNodes,9);  % right hand side (3 stress components + 3 strain components + 3 moments)
M = sparse(numNodes,numNodes);  

[qp,w] = Gauss(nnpe);
for e = 1:numEls
    % Elemental initialization
    Me = zeros(nnpe);        
    Re = zeros(nnpe,9);
    
    % get nodal coordinates for this element
    xvec=nCoords(elCon(e,:),1); % these are column vectors of the
    yvec=nCoords(elCon(e,:),2); % physical coordinates of the
                                    % nodes of this element

    x_elemlength =  abs(xvec(2) - xvec(1))/2;                           
    y_elemlength =  abs(yvec(3) - yvec(2))/2;         
    
    
    lcU = glU(gatherMat(e,:)); % local displacement
    for iqp = 1:length(w)
        % Calcluate N, dN, ddN matrices
        NofXiEta=StrainStressInterp(qp(iqp,:));    % 1x4 interpolation matrix
        ddNdXiEta=ddNmatrix(qp(iqp,:)); % 3x12 matrix, of second derivs of shape func
                                        % at this QP in parent domain
        
        % jacobian at this QP (2x2 matrix)
        JofXiEta= [x_elemlength, 0; 0, y_elemlength];
        
        % now create the B matrix
        B =  [(1/x_elemlength)^2*ddNdXiEta(1,:); (1/y_elemlength)^2*ddNdXiEta(2,:); ...
            (2/x_elemlength)*(1/y_elemlength)*ddNdXiEta(3,:)];
        
        % elemental strain and stress
        lcStrain = - layer_position * B * lcU;
        lcStress = D * lcStrain;
        lcMoment = - D * thickness^3 * B * lcU / 12;
        
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

strain = M\R(:,1:3);
stress = M\R(:,4:6);
moment = M\R(:,7:9);
end

function N=StrainStressInterp(LocPos)
% N=Nmatrix(QP)
% This returns the matrix of shape functions N calculated at the local
% coordinates XI,ETA for 2D elements.
% last edit: 29 April 2015 H. Ritz

xi=LocPos(1); eta=LocPos(2);
N = 1/4*[(1 - xi).*(1-eta),(1+xi).*(1-eta),...
          (1+xi).*(1+eta), (1-xi).*(1+eta)];  % Calculate the 4 basis functions 
                                              % at (xi,eta). N is a row
                                              % vector
end


