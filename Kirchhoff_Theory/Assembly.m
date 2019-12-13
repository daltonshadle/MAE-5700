% *************************************************************************
% Name: Assembly.m
% Notes: globalSystem=Assembly(MeshStruct);
%     Assemble global stiffness matrix K and global force vector F for the
%     LinElast code. 
%    last edit: 22 November 2017 H. Ritz
% Project Updates: Pulled from LinElast code base from MAE-5700 course,
%     updates were made to implement Kirchhoff Theory Plate Bending
% Update Authors: Dalton and Sairam
% *************************************************************************

function globalSystem=Assembly(MeshStruct)
%% unpack necessary information
nnpe=MeshStruct.nnpe;
numDOF=MeshStruct.numDOF;
numEls=MeshStruct.numEls;
numEq=MeshStruct.numEq;
gatherMat=MeshStruct.gatherMat;

%% initialize the global system
K=zeros(numEq,numEq); 
F=zeros(numEq,1);


% get the appropriate quadrature locations and weights
[qp,w] = Gauss(nnpe);  % quadrature points

% iterate over each element to assemle elemental force and stiffness
% matrix, add to global force and stiffness matrix
for e=1:numEls
    % make the local stiffness matrix and local force vector for this element
    [ke,fe] = TwoDElem(MeshStruct, e, qp, w); 
    
    for Lrw = 1 : (nnpe*numDOF)
        % global row index
        Grw = gatherMat(e,Lrw);
        
        % add elemetal force to global force
        F(Grw) = F(Grw) + fe(Lrw);
        for Lcl = 1 : (nnpe*numDOF)
            % global column index
            Gcl = gatherMat(e,Lcl);
            
            % Assemble local stiffness matrix into the global
            % stiffness matrix here
            K(Grw,Gcl) = K(Grw,Gcl) + ke(Lrw,Lcl);
        end
    end
end

%% Package variables into the output structs. 
% NOTE: both K and F are made sparse to greatly reduce solution time.
globalSystem.K=sparse(K);
globalSystem.F=sparse(F);
