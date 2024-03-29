function globalSystem=Assembly(globalSystem,meshStruct)
% globalSystem=ASSEMBLY(globalSystem,meshStruct)
% Assemble global stiffness matrix K for the TRUSS2D code. 
% last edit: 5 September 2017 H. Ritz

% unpack necessary inputs
nnpe=meshStruct.nnpe;
numDOF=meshStruct.numDOF;
numEls=meshStruct.numEls;
gatherMat=meshStruct.gatherMat;

% We will first create a local stiffness matrix "ke" (call the function 
% TrussElem.m) and then put "ke" into our global stiffness matrix K

% The following nested for loop may seem intimidating, but all we need to
% do is to locate the local and global row & column index.
numNodes=meshStruct.numNodes;
elCon=meshStruct.elCon;
K=zeros(numNodes*numDOF);
for e=1:numEls % for each element
    ke = TrussElem(e,meshStruct); % make the local stiffness matrix
    for i = 1 : nnpe
        for ii = 1 : numDOF
            Lrw = numDOF*(i-1)+ii; % local row index
            Grw = numDOF*(elCon(e,i)-1)+ii; % global row index
            for j = 1 : nnpe
                for jj = 1 : numDOF
                    Lcl = numDOF*(j-1)+jj; % local column index
                    Gcl = numDOF*(elCon(e,j)-1)+jj; % global column index
                    
                    % Assemble local stiffness matrix into the global
                    % stiffness matrix here
                    K(Grw,Gcl) = K(Grw,Gcl) + ke(Lrw,Lcl); 
                    
                end
            end
        end
    end
end

% Package variables into the output structs
globalSystem.K=K;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comment out return to enact assignment errors.
%return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here is another version of the same code. You can take advantage of the
% gatherMat array which has already been defined to simplify some of the
% loops. Fill in the missing code to assemble ke into K2 which should be
% identical to K. 

K2 = zeros(numNodes*numDOF);
for e = 1:numEls % for each element
    ke = TrussElem(e, meshStruct); % make the local stiffness matrix
    for lclROW = 1:(nnpe*numDOF)
        glbROW = gatherMat(e, lclROW); % global row index
        % need to find the lclROW number in terms of the glbROW number from
        % the assembly matrix, or in this case gatherMat
        
        for lclCOL = 1:(nnpe*numDOF)
            glbCOL = gatherMat(e, lclCOL); % global column index
            % need to find the lclCOL number in terms of the glbCOL number 
            % from the assembly matrix, or in this case gatherMat
            
            % Assemble local stiffness matrix into the global
            % stiffness matrix here
            K2(glbROW,glbCOL) = K2(glbROW,glbCOL) + ke(lclROW, lclCOL);
            % need to assign K2 at the global node numbers to ke at the 
            % local node numbers
        end
    end

end

% Compare results from different methods.
display(['The maximum difference between K and K2 is ', num2str(max(max(abs(K-K2))))]);

% *********************** Multipoint Constraint Code **********************
% Pull all variables from globalSystem struct
c = globalSystem.c;
multiConstraint = globalSystem.multiConstraint;

% Initialize numMPC to the number of multipoint constraint
numMPC = size(multiConstraint, 1);
globMPC = [];

% iterate over each multipoint constraint and construct c matrix
for i = 1:numMPC
    % intialize two MPC variables, node number and incline angle
    nodeNumMPC = multiConstraint(i, 1);
    angleMPC = multiConstraint(i, 2);
    
    % create tempC matrix holding the local node MPC values
    tempC = [sin(angleMPC), -cos(angleMPC)];
    
    % add to correct location in global c
    globMPCmap = [nnpe * (nodeNumMPC-1) + 1, nnpe * (nodeNumMPC-1) + 2];
    globMPC = [globMPC, globMPCmap];
    c(i, globMPCmap) = c(i, globMPCmap) + tempC;
end
globMPC = unique(globMPC);

% Reassign c to globalSystem struct to pass back through
globalSystem.c = c;
globalSystem.globMPC = globMPC;



