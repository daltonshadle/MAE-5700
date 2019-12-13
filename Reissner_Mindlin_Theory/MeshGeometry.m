% *************************************************************************
% Name: MeshGeometry.m
% Notes: [meshStruct,boundStruct,PlotInstructions]=MeshGeometry;
%     Define the mesh and identify the boundary nodes and elements. See the
%     help under LINELAST for a description of the output structs. Note: the
%     only change from TwoDBVP version of this function is numDOF=2;
%     last update: 15 Nov 2015 H. Ritz; Y. Xu
% Project Updates: Pulled from LinElast code base from MAE-5700 course,
%     updates were made to implement Kirchhoff Theory Plate Bending
% Update Authors: Dalton and Sairam
% *************************************************************************

function [meshStruct,boundStruct,PlotInstructions]=MeshGeometry
%% Plot instruction
PlotInstructions.plot_mesh     = 'yes';
PlotInstructions.plot_node     = 'yes';
PlotInstructions.plot_boundary = 'yes';
PlotInstructions.plot_contour  = 'yes';
PlotInstructions.plot_vector   = 'yes';
PlotInstructions.plot_deformed = 'yes';
PlotInstructions.plot_strains  = 'yes';
PlotInstructions.plot_fringes  = 'yes';

%% Input parameters
% number of nodes per element, currently Q4 elements are supported
nnpe = 4;        

% number of spatial dimensions
nsd = 2;

% Domain
xl = 0.0;           % left location of the range in the x direction
xr = 1.0;           % right location of the range in the x direction
yb = 0.0;           % bottom location of the range in the y direction
yt = 1.0;           % top location of the range in the y direction

% thickness
thickness = (xr-xl)/10;
zb = -thickness/2;  % bottom location of the range in the z direction
zt = thickness/2;   % top location of the range in the z direction
thickness_layer = linspace(zb, zt, 30);

% number of elements in x and y dirns
nx = 40; ny = 40;

% Use the provided box grid generator
[nCoords,elCon,boundStruct]=BoxGrid_2D(nsd,xl,xr,yb,yt,nnpe,nx,ny);

%% Connectivity array
% Use the connectivity array to define gatherMat which shows the global
% degrees of freedom for each local degree of freedom.
numEls = size(elCon,1);
numNodes=size(nCoords,1);
numDOF = 3; % linear elasticity is a vector field problem
gatherMat=zeros(numEls,(nnpe*numDOF));
for n=1:nnpe  % loop over the number of nodes per element
    globalNodes=elCon(:,n); % global node number for this local node
    for d=1:numDOF % loop over the number of degrees of freedom per node
        % use global node numbers to find global DOFs
        % corresponding to the local DOFs.
        gatherMat(:,(n-1)*numDOF+d)=(globalNodes-1)*numDOF+d;
    end
end

%% Package the necessary variables into the output structure
meshStruct.nnpe       = nnpe;
meshStruct.nCoords    = nCoords;
meshStruct.elCon      = elCon;
meshStruct.nsd        = nsd;
meshStruct.numNodes   = numNodes;
meshStruct.numEls     = numEls;
meshStruct.numDOF     = numDOF;
meshStruct.numEq      = numDOF*numNodes;
meshStruct.gatherMat  = gatherMat;
meshStruct.thickness  = thickness;
meshStruct.thickness_layer = thickness_layer;
meshStruct.nx         = nx;
meshStruct.ny         = ny;

%% Plot the mesh
% PlotGrid(meshStruct,boundStruct,PlotInstructions);
