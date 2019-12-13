function [meshStruct,boundStruct,PlotInstructions]=MeshGeometry(nx,ny,nnpe)
% [meshStruct,boundStruct,PlotInstructions]=MeshGeometry; 
% Define the mesh and identify the boundary nodes and elements. See the
% help under LINELAST for a description of the output structs. Note: the
% only change from TwoDBVP version of this function is numDOF=2;
%
% last update: 15 Nov 2015 H. Ritz; Y. Xu


% Choose a meshing method
% "1" for mesh directly generated from Matlab "BoxGrid_2D" function
% "2" for mesh generated from Ansys and loaded from "loadFromGridFile"
% function
LoadChoice = 1;

% Plot instruction
PlotInstructions.plot_mesh     = 'no';  % What to plot. For big meshes, 
PlotInstructions.plot_node     = 'no';  % it is better not to plot node and vector.
PlotInstructions.plot_boundary = 'no';  % Change this information as appropriate
PlotInstructions.plot_contour  = 'no';
PlotInstructions.plot_vector   = 'no';
PlotInstructions.plot_deformed = 'no';
PlotInstructions.plot_fringes  = 'no';

nnpe = 3;           % number of nodes per element.
                    % currently T3 and Q4 elements are supported

nsd = 2;            % number of spatial dimensions 
if LoadChoice==1 % use LoadChoice =1 only when the geometry is a rectangular
    % Set this information for each problem

    xl = 0.0;           % left location of the range in the x direction
    xr = 1.0;           % right location of the range in the x direction
    yb = 0.0;           % bottom location of the range in the y direction
    yt = 0.01;           % top location of the range in the y direction

    nx = 79*3;            % number of elements in x direction
    ny = 1;            % number of elements in y direction

    % Use the provided box grid generator
    [nCoords,elCon,boundStruct]=BoxGrid_2D(nsd,xl,xr,yb,yt,nnpe,nx,ny);
    
    % Mesh Generator for HW12 Problem 1 (NOTE: nx must be multiple of 3)
    % [nCoords,elCon,boundStruct]=BoxGrid_2D_P1(nsd,xl,xr,yb,yt,nnpe,nx,ny);
elseif LoadChoice==2
    [nCoords,elCon,boundStruct]=loadFromGridFile('ansys');
else
    error('Invalid LoadChoice');
end

% Use the connectivity array to define gatherMat which shows the global
% degrees of freedom for each local degree of freedom.

numEls = size(elCon,1);
numNodes=size(nCoords,1);
numDOF = 2; % linear elasticity is a vector field problem
gatherMat=zeros(numEls,(nnpe*numDOF));
for n=1:nnpe  % loop over the number of nodes per element
    globalNodes=elCon(:,n); % global node number for this local node
    for d=1:numDOF % loop over the number of degrees of freedom per node
        % use global node numbers to find global DOFs
        % corresponding to the local DOFs.
        gatherMat(:,(n-1)*numDOF+d)=(globalNodes-1)*numDOF+d;
    end
end

% package the necessary variables into the output structure
meshStruct.nnpe     = nnpe;
meshStruct.nCoords  = nCoords;
meshStruct.elCon    = elCon;
meshStruct.nsd      = nsd;
meshStruct.numNodes = numNodes;
meshStruct.numEls   = numEls;
meshStruct.numDOF   = numDOF; 
meshStruct.numEq    = numDOF*numNodes;
meshStruct.gatherMat= gatherMat;
meshStruct.yt = yt;

% Plot the mesh
PlotGrid(meshStruct,boundStruct,PlotInstructions);
