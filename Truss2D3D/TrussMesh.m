function meshStruct=TrussMesh;
% meshStruct=TRUSSMESH
% This file defines the geometry for the TRUSS2D3D code. User must input the
% nodal coordinates and the connectivity array.
% The portions of the code that might change for each new problem are
% clearly indicated. The output from this function is meshStruct which has
% all of the mesh information, as detailed in the help documentation for
% Truss2D.

% last edit: 30 July 2015 H. Ritz

% some information about the types of elements
nnpe=2;   % number of nodes per element. define this variable since it will 
          % be used for determining the dimensions of future arrays.
        
% Define the coordinates of the nodes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE THIS FOR EACH PROBLEM
% If you give an x and y coordinate, this is a 2D problem. If you give 3
% coordinates for each node, this is a 3D problem.

% first column is x coordinates of each node
% second column is y coordinates of each node
% third column (if present) is z coordinates of each node
nCoords = [2 0 0;
           0 0 1;
           0 0 3;
           0 2 1;];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
numNodes=size(nCoords,1); % number of nodal points in the mesh. 
                          % does not have to be entered manually.
numDim=size(nCoords,2); % define this variable to be able to 
                        % change code easily for 2D or 3D problems.
numDOF=numDim; % number of degrees of freedom per node. 

% Define the connectivity array (dimensions are numEls X nnpe)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE THIS FOR EACH PROBLEM

% elCon(i,j) is the global node number of the 
% jth node of the ith element
elCon = [1 2;
         1 3;
         1 4;];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numEls=size(elCon,1); % number of elements in the mesh.
                      % does not have to be entered manually.

% Use the connectivity array to define gatherMat which shows the global
% degrees of freedom for each local degree of freedom.
gatherMat=zeros(numEls,(nnpe*numDOF));
for n=1:nnpe  % loop over the number of nodes per element
    globalNodes=elCon(:,n); % global node number for this local node
    for d=1:numDOF % loop over the number of degrees of freedom per node
        % use global node numbers to find global DOFs
        % corresponding to the local DOFs.
        gatherMat(:,(n-1)*numDOF+d)=(globalNodes-1)*numDOF+d;
    end
end

% Package variables into the mesh struct
meshStruct.nCoords  =nCoords; % global nodal coordinates
meshStruct.elCon    =elCon;   % connectivity array
meshStruct.nnpe     =nnpe;    % number of nodes per element (2)
meshStruct.numDim   =numDim;  % number of spatial dimensions (2 or 3)
meshStruct.numDOF   =numDOF;  % degrees of freedom per nodes
meshStruct.numNodes =numNodes;% number of nodes
meshStruct.numEls   =numEls;  % number of elements
meshStruct.gatherMat=gatherMat;% global DOF for local DOF
