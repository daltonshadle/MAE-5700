function meshStruct=TrussMesh;
% meshStruct=TRUSSMESH
% This file defines the geometry for the TRUSS2D code. User must input the
% nodal coordinates and the connectivity array.
% The portions of the code that might change for each new problem are
% clearly indicated. The output from this function is meshStruct which has
% all of the mesh information, as detailed in the help documentation for
% Truss2D.

% last edit: 10 July 2015 H. Ritz

% some information about the types of elements
nnpe=2;   % number of nodes per element. define this variable since it will 
          % be used for determining the dimensions of future arrays.
        
% Define the coordinates of the nodes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE THIS FOR EACH PROBLEM
nCoords=[  0   0; % first column is x coordinates of each node
         0.5   0; % second column is y coordinates of each node
         0.5 0.5];
% Variable for question 6
nCoords=[  0,   0;
         0.5,   0;
           1,   0;
         1.5,   0;
           2,   0;
           0, 0.5;
         0.5, 0.5;
           1, 0.5;
         1.5, 0.5;
           2, 0.5];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the connectivity array (dimensions are numEls X nnpe)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE THIS FOR EACH PROBLEM
elCon=[1 3;  % elCon(i,j) is the global node number of the 
       2 3;  % jth node of the ith element
       1 2];

% Variable for question 6
elCon=[1, 6;
       1, 7;
       2, 7;
       2, 8;
       3, 8;
       4, 8;
       4, 9;
       5, 9;
       5, 10;
       1, 2;
       2, 3;
       3, 4;
       4, 5;
       6, 7;
       7, 8;
       8, 9;
       9, 10];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
numNodes=size(nCoords,1); % number of nodal points in the mesh. 
                          % does not have to be entered manually.
numDim=size(nCoords,2);   % define this variable to be able to 
                          % change code easily for 2D or 3D problems.
numDOF=numDim;            % number of degrees of freedom per node. 


numEls=size(elCon,1);     % number of elements in the mesh.
                          % does not have to be entered manually.

% Use the connectivity array to define gatherMat
% gatherMat(i,j) is the global degree of freedom for element i's jth local
% degree of freedom
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
