function meshStruct=BarMesh;
% meshStruct=BARMESH 
% This file defines the geometry for ONEDBVP finite element code. You must
% choose the element type, number of elements, and number of quadrature
% points for evaluating integrals. Also set the physical domain for the
% problem.
%
% This function calculates the nodal coordinates, element connectivity
% array, and gatherMat between local and global degrees of freedom. 


% last edit: 5 August 2015 H. Ritz

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define this for each problem
nnpe=3;   % number of nodes per element: 
          % 2 for linear elements
          % 3 for quadratic elements
          % 4 for cubic elements
domain=[0 1]; % the x values of the limits of the domain
% domain=[2 6]; % the x values of the limits of the domain
numEls=50; % number of elements in the mesh
numQP=2;    % number of quad points for numerical integration. 
            % choose this to be sufficient for whatever integrals 
            % you need to evaluate.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numDOF=1; % number of degrees of freedom per node.
          % some routines need this variable defined.

% use the domain, number of elements, and element type to find the nodal
% point locations and make the connectivity array.
xMin=domain(1); xMax=domain(2);
switch nnpe % depending on element type ...
    case 2 % 2 nodes per element, or linear elements
        numNodes=numEls+1; % for linear elements there is one more node 
                           % than the number of elements
        nCoords=linspace(xMin,xMax,numNodes); % coordinates of the 
                                              % evenly spaced nodes
        elCon=[1:numNodes-1;2:numNodes]'; % element connectivity array
    case 3
        numNodes=2*numEls+1; % for quad elements
        nCoords=linspace(xMin,xMax,numNodes); % coordinates of the 
                                              % evenly spaced nodes
        elCon=[1:2:numNodes-2; 2:2:numNodes-1; 3:2:numNodes]'; % element connectivity array
    case 4
        numNodes=3*numEls+1; % for quad elements
        nCoords=linspace(xMin,xMax,numNodes); % coordinates of the 
                                              % evenly spaced nodes
        elCon=[1:3:numNodes-3; 2:3:numNodes-2; 3:3:numNodes-1; 4:3:numNodes]'; % element connectivity array
    otherwise
        error('Invalid element type.')
end

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
numEq=numNodes*numDOF; % number of equations in the global system

% Package variables into the mesh struct
meshStruct.nCoords  =nCoords;  % global nodal coordinates
meshStruct.elCon    =elCon;    % connectivity array
meshStruct.nnpe     =nnpe;     % number of nodes per element (2)
meshStruct.numDOF   =numDOF;   % degrees of freedom per node
meshStruct.numNodes =numNodes; % number of nodes
meshStruct.numEls   =numEls;   % number of elements
meshStruct.numQP    =numQP;    % number of quad points for integrals
meshStruct.numEq    =numEq;    % number of equations in global system
meshStruct.gatherMat=gatherMat;% global DOF for local DOF



