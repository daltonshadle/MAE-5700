function nodalGradient = calNodalGradient(glField,meshStruct)
%nodalGradient = calNodalGradient(glField,meshStruct)   Calculate nodal
%gradient values 
%
% Input 
% glField    : column or row vector, the global nodal solution of 2DBVP;
% meshStruct : structure containing FE mesh info, see TwoDBVP.m for
% detailed spec.
%
% Output
% nodalGradient : nx2 matrix, n is the number of node. Nodal gradient
% values. First column dT/dx, second column dT/dy.

% Algorithm
% This code first calculates gradient values at gauss points, and
% then extrapolates to the nodal points. That way the output gradients
% are one order more accurate than calculating them directly at nodal
% points. This technique in FE is called "stress recovery", and can be
% applied in both scalar field (2DBVP) and vector field (Linear elastic)
% problems.
% The idea of the extrapolation is global least square. In the end, we form
% M*u=R, where M is the global mass matrix M = integral(N'*N) and 
% R = integral(N'*filed_gp). u is the nodal gradient values. See lecture
% notes for more info.

% last edit: Nov 9, 2016 Y. Xu

if isrow(glField)
    glField = glField'; % to column vector
end

% unpack what you need
numNodes = meshStruct.numNodes;
nnpe=meshStruct.nnpe;
numEls = meshStruct.numEls;
nCoords=meshStruct.nCoords;
elCon=meshStruct.elCon;
gatherMat=meshStruct.gatherMat;

% Initialization
R = zeros(numNodes,2); % right hand side (first column dT/dx, 
                       % second column dT/dy)
M = spalloc(numNodes,numNodes,10*numNodes); % allocate sparse matrix, the 
                                            % last input is the estimation
                                            % of nonzero numbers of M
[qp,w] = Gauss(nnpe);
for e=1:numEls

    xynode=nCoords(elCon(e,:),:); % nodal coordinate
    lcField=glField(gatherMat(e,:)); % from global to local nodal solutions
    
    % Elemental initialization 
    Me = zeros(nnpe);
    Re = zeros(nnpe,2);
    for iqp=1:length(w) % loop over quadrature points
        
        NofXiEta=Nmatrix(qp(iqp,:),nnpe); % row vector of shape functions
        % at this QP in parent domain
        dNdXiEta=dNmatrix(qp(iqp,:),nnpe); % row vector of shape function derivatives
        % at this QP in parent domain
        JofXiEta=dNdXiEta*xynode; % jacobian at this QP 
        B=JofXiEta\dNdXiEta;% row vector of dNdX at ths QP, needed for integration
        
        % elemental gradient at quad point
        gaussGradientEle=B*lcField; 
        
        % M matrix and right hand side for global least square recovery of
        % the gradient field
        Me = Me + NofXiEta'*NofXiEta * w(iqp) * det(JofXiEta);
        Re = Re + NofXiEta'*gaussGradientEle' * w(iqp) * det(JofXiEta);

    end
    
    % Assembly: we can do the same thing as in Assembly.m, but here we
    % use a simple way--a little inefficient for sparse matrix. This would
    % still be better than (when solving the large linear system) using
    % full matrix.
    
    % global index
    glInd = gatherMat(e,:);
    
    % Assemble M and R
    M(glInd,glInd) = M(glInd,glInd) + Me;
    R(glInd,:)  = R(glInd,:) + Re;
end

nodalGradient = M\R;
