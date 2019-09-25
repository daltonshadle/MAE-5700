function globalSystem=Assembly(globalSystem,meshStruct)
% globalSystem=ASSEMBLY(globalSystem,meshStruct)
% Assemble global stiffness matrix K for the TRUSS2D3D code. 
% last edit: 30 July 2015 H. Ritz

% unpack necessary inputs
nnpe=meshStruct.nnpe;
numDOF=meshStruct.numDOF;
numNodes=meshStruct.numNodes;
numEls=meshStruct.numEls;
gatherMat=meshStruct.gatherMat;
K=globalSystem.K;

% for each element, make the local stiffness matrix and then assemble
% "ke" into our global stiffness matrix K
f_total_therm = zeros(numNodes * numDOF, 1);
for e=1:numEls % for each element
    [ke, fe_therm] = TrussElem(e,meshStruct); % make the local stiffness matrix
    for Lrw = 1 : (nnpe*numDOF)
        Grw = gatherMat(e,Lrw); % global row index
        for Lcl = 1 : (nnpe*numDOF)
            Gcl = gatherMat(e,Lcl); % global column index
            
            % Assemble local stiffness matrix into the global
            % stiffness matrix here
            K(Grw,Gcl) = K(Grw,Gcl) + ke(Lrw,Lcl);
            
        end
    end
	local2global_map = gatherMat(e, :);
	f_total_therm(local2global_map) = f_total_therm(local2global_map) + fe_therm;	
end

% Package variables into the output structs
globalSystem.K=K;
globalSystem.f_total_therm = f_total_therm;
