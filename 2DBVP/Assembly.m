function globalSystem=Assembly(MeshStruct)
% globalSystem=ASSEMBLY(meshStruct)
% Assemble global stiffness matrix K and global force vector F for the
% TWODBVP code. 

% last edit: 7 November 2017 H. Ritz

% unpack necessary information
nnpe  =MeshStruct.nnpe;
numDOF=MeshStruct.numDOF;
numEls=MeshStruct.numEls;
numEq =MeshStruct.numEq;
gatherMat=MeshStruct.gatherMat;

% initialize the global system
K=zeros(numEq,numEq);
F=zeros(numEq,1);

% for each element, make the local stiffness matrix and local force vector
% then assemble "ke" into our global stiffness matrix K and "fe" into our
% global force vector F

% get the appropriate quadrature locations and weights
[qp,w] = Gauss(nnpe);
for e=1:numEls
    [ke,fe] = TwoDElem(MeshStruct, e, qp, w); % make the local stiffness matrix and 
                                              % local force vector for this element
    
    for Lrw = 1 : (nnpe*numDOF)
        Grw = gatherMat(e,Lrw); % global row index
        F(Grw) = F(Grw) + fe(Lrw);
        for Lcl = 1 : (nnpe*numDOF)
            Gcl = gatherMat(e,Lcl); % global column index
            
            % Assemble local stiffness matrix into the global
            % stiffness matrix here
            K(Grw,Gcl) = K(Grw,Gcl) + ke(Lrw,Lcl);
            
        end
    end
end

% Package variables into the output structs. NOTE: both K and F are made
% sparse to greatly reduce solution time.
globalSystem.K=sparse(K);
globalSystem.F=sparse(F);
