function globalSystem=Assembly(globalSystem,meshStruct)
% globalSystem=ASSEMBLY(globalSystem,meshStruct)
% Assemble global stiffness matrix K and global force vector F for the
% FRAME code. This function makes use of the MATLAB command SPARSE which
% takes as input three vectors and assembles a sparse matrix from the
% information contained.
% last edit: 27 September 2017 H. Ritz

% unpack necessary inputs
nnpe=meshStruct.nnpe;
numDOF=meshStruct.numDOF;
numEls=meshStruct.numEls;
gatherMat=meshStruct.gatherMat;
K=globalSystem.K;
F=globalSystem.F;

% for each element, make the local stiffness matrix and local force vector
% then assemble "ke" into our global stiffness matrix K and "fe" into our
% global force vector
for e=1:numEls
    [ke,fe] = FrameElem(e,meshStruct); % make the local stiffness matrix 
                                       % and local force vector for 
                                       % this element
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

