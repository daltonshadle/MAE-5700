% *************************************************************************
% Name: ApplyNaturalBC.m
% Notes: globalSystem=ApplyNaturalBC(i,boundStruct,meshStruct,globalSystem)
%     This function applies natural BCs
%     N. Zabaras (zabaras@cornell.edu) & Xiang Ma (xm25@cornell.edu) |
%     Last update: 17 Nov 2015 H. Ritz ; Y. Xu
% Project Updates: Pulled from LinElast code base from MAE-5700 course,
%     updates were made to implement Kirchhoff Theory Plate Bending
% Update Authors: Dalton and Sairam
% *************************************************************************

function  globalSystem=ApplyNaturalBC(i,boundStruct,meshStruct,globalSystem)
%% unpack necessary variables
sideInd = boundStruct.SurfNat(i,1);
fDOF = boundStruct.SurfNat(i,2); % force DOF
q = boundStruct.SurfNat(i,3); % force value
BElems = boundStruct.elements(sideInd).Elems;
SurfID = boundStruct.elements(sideInd).SurfaceIndicator;

nnpe=meshStruct.nnpe;
numDOF=meshStruct.numDOF;
nCoords=meshStruct.nCoords;
elCon=meshStruct.elCon;
gatherMat=meshStruct.gatherMat;
F = globalSystem.F;
ldof = nnpe*numDOF;       % Degrees of freedom per element

Num = length(BElems); %  Extract the number of boundary elements

% get the appropriate quadrature locations and weights
[qp,w] = OneDGauss(nnpe);  % quadrature points

% construct q vector
q_vec = zeros(3,1);
q_vec(fDOF,1) = q;

% iterate over each element
for elmID = 1 : Num
    fe = zeros(ldof,1);   % initialize element load vector
    glb = elCon(BElems(elmID),:);   % Global number of the element nodes
    
    % get nodal coordinates for this element
    xvec=nCoords(glb,1); % these are column vectors of the
    yvec=nCoords(glb,2); % physical coordinates of the nodes of this element
    coord = [xvec,yvec];
    
    x_elemlength =  abs(xvec(2) - xvec(1))/2;
    y_elemlength =  abs(yvec(3) - yvec(2))/2;
    
    % iterate over each quadrature point
    for iqp=1:size(qp,1)
        % Calcluate N, dN, ddN matrices
        switch (SurfID(elmID))
            case -2
                NofXiEta=Nmatrix([qp(iqp), -1]); % 1x12 matrix,  row vector of shape functions
                                                 % at this QP in parent domain
                dNdXiEta=dNmatrix([qp(iqp), -1]);  % 2x12 matrix, of first derivs of shape func
                                                   % at this QP in parent domain
                % jacobian at this QP (1x2 matrix)
                JofXiEta= [x_elemlength, 0];
            case +1
                NofXiEta=Nmatrix([1, qp(iqp)]); % 1x12 matrix,  row vector of shape functions
                                                 % at this QP in parent domain
                dNdXiEta=dNmatrix([1, qp(iqp)]);  % 2x12 matrix, of first derivs of shape func
                                                   % at this QP in parent domain
                
                % jacobian at this QP (1x2 matrix)
                JofXiEta= [0, y_elemlength];
            case +2
                NofXiEta=Nmatrix([qp(iqp), 1]); % 1x12 matrix,  row vector of shape functions
                                                 % at this QP in parent domain
                dNdXiEta=dNmatrix([qp(iqp), 1]);  % 2x12 matrix, of first derivs of shape func
                                                   % at this QP in parent domain
                % jacobian at this QP (1x2 matrix)
                JofXiEta= [x_elemlength, 0];
            case -1
                NofXiEta=Nmatrix([-1, qp(iqp)]); % 1x12 matrix,  row vector of shape functions
                                                 % at this QP in parent domain
                dNdXiEta=dNmatrix([-1, qp(iqp)]);  % 2x12 matrix, of first derivs of shape func
                                                   % at this QP in parent domain
                
                % jacobian at this QP (1x2 matrix)
                JofXiEta= [0, y_elemlength];
        end
        
        Jacc = dNdXiEta*coord;      % Calculate the components of the Jaccobian matrix
        
        % now create the N 
        N=[];
        for np=1:nnpe
            N=[N,[NofXiEta(np), 0, 0; 0, NofXiEta(np), 0; 0, 0, NofXiEta(np)]];
        end
    
        % add the weighted contributions of this QP to the elemental stiffness
        % matrix and elemental body force vector
        fe=fe+(N'*q_vec)*w(iqp)*norm(JofXiEta);
    end
    
    % Assemble elemental force vector to global force vector
    % global row indices
    glbROW = gatherMat(BElems(elmID),:);
    F(glbROW)=F(glbROW)+fe;
end

%% Package variables in structs
globalSystem.F = F;