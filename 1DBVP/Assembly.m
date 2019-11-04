function globalSystem=Assembly(meshStruct, boundStruct)
% globalSystem=ASSEMBLY(meshStruct, boundStruct)
% Assemble global stiffness matrix K and global force vector F for the
% ONEDBVP code. 

% last edit: 5 August 2015 H. Ritz

% unpack necessary inputs
nnpe  =meshStruct.nnpe;
numDOF=meshStruct.numDOF;
numEq =meshStruct.numEq;
numQP =meshStruct.numQP;
numEls=meshStruct.numEls;
gatherMat=meshStruct.gatherMat;

% initialize the global system
K=zeros(numEq);
F=zeros(numEq,1);


% for each element, make the local stiffness matrix and local force vector
% then assemble "ke" into our global stiffness matrix K and "fe" into our
% global force vector F

% get the appropriate quadrature locations and weights
[qp,w] = Gauss(numQP);
for e=1:numEls
    % make the local stiffness matrix and local force vector for this element
    [ke,fe] = OneDElem(e, qp, w, meshStruct, boundStruct);   
    
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

% Package variables into the output structs
globalSystem.K=K;
globalSystem.F=F;
