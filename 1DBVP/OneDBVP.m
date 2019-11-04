% OneDBVP
% Finite Element code for solving 1D boundary value problems, such as
% linear elasticity or 1D heat conduction. The code will solve any problem
% that can be written as:
% -d/dx( p(x)du/dx )=f(x)
% 
% This program supports linear, quadratic, and cubic elements, as well as
% using gauss quadrature to evaluate all integrals in the weak form.
% To change the problem being analyzed, edit the following functions:
% -BarMesh
% -InputData
% -FF, PP
% -Exact
% All variables are stored in structures, or "structs" which are detailed
% here:
% meshStruct.nCoords   % global nodal coordinates
% meshStruct.elCon     % connectivity array
% meshStruct.nnpe      % number of nodes per element (2)
% meshStruct.numDOF    % degrees of freedom per node
% meshStruct.numNodes  % number of nodes
% meshStruct.numEq     % number of equations in global system
% meshStruct.numEls    % number of elements
% meshStruct.gatherMat % global DOF for local DOF
% 
% globalSystem.K           % global stiffness matrix
% globalSystem.F           % global force vector including reaction vector
% globalSystem.d           % global solution vector including essential BCs
% globalSystem.reactionVec % reaction vector on essential DOF
% globalSystem.displacement% transverse displacement at many points along
                           % the beam
% 
% boundStruct.boundCond  % information about BCs
% boundStruct.appForces  % information about applied forces


% last edit: 23 October 2017 H. Ritz

% Preliminary steps
home; clear; close all;  % clean the workspace

% Preprocessing
meshStruct=BarMesh;       % make the mesh geometry
boundStruct=InputData(meshStruct); % read problem data about 
                                   % applied forces and BCs 

% Calculation and assembly of global stiffness matrix
globalSystem=Assembly(meshStruct, boundStruct);  % creates local stiffness  
                                    % matrices and uses connectivity array 
                                    % to create global stiffness matrix. 

% Solution Phase
[globalSystem,boundStruct]=...
    ApplyBC(globalSystem, meshStruct, boundStruct); % applies both natural 
                                                    % and essential BCs
globalSystem=...
    Soln(globalSystem,meshStruct,boundStruct); % solves the global system 
                                               % for nodal DOF 
                                               % and reaction "forces"

% Postprocessor Phase 
globalSystem=...
    PostProcess(globalSystem,meshStruct); % plot u(x), du/dx, 
                                          % analyze error, etc.

                                                      
                                                      