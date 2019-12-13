% *************************************************************************
% Name: HW12_Q1.m
% Authors: Dalton and Sairam
% Date: 12/7/19
% Notes: Summarized code changes for HW12 Question 1
% *************************************************************************

% From BoxGrid_2D.m *******************************************************
% ...
% Add extra boundary surface for bottom 2 boundary conditions, Surface 1 is
% our traction free surface and Surface 5 is our symmetry condition
BoundaryNodes(1).Nodes = [1:(nelemx)/3];  % The first boundary (bottom left 1/3)
BoundaryNodes(5).Nodes = [(nelemx)/3+1:nelemx+1];  % The fifth boundary (bottom right 2/3)
%...
% Add extra boundary surface for bottom 2 boundary conditions, Surface 1 is
% our traction free surface and Surface 5 is our symmetry condition
BoundaryElems(1).Elems = [1:nelemx/3];  % The first boundary (bottom left 1/3)
BoundaryElems(1).SurfaceIndicator = -2*ones(nelemx/3,1);

BoundaryElems(5).Elems = [nelemx/3+1:nelemx];  % The fifth boundary (bottom right 2/3)
BoundaryElems(5).SurfaceIndicator = -2*ones(2*nelemx/3,1);


% From InputData.m *******************************************************
% ...
% Essential BC's for HW12 Problem 1, symmetry displacement BC's on surface 4,5
boundStruct.SurfEssV = [5 2 0;
                        4 1 0];
% ...
% Natural BC's for HW12 Problem 1, traction free surfaces 1,2, normal
% traction on surface 3
boundStruct.SurfNat = [3 0 50e6;
                       2 0 0;
                       1 0 0];

% From MeshGeomtry.m ******************************************************
% ...
xl = 0.0;           % left location of the range in the x direction
xr = 0.15;           % right location of the range in the x direction
yb = 0.0;           % bottom location of the range in the y direction
yt = 0.5;           % top location of the range in the y direction

% Use the provided box grid generator
% Mesh Generator for HW12 Problem 1 (NOTE: nx must be multiple of 3)
[nCoords,elCon,boundStruct]=BoxGrid_2D(nsd,xl,xr,yb,yt,nnpe,nx,ny);













