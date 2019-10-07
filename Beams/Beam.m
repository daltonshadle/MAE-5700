% BEAM
% Finite Element code for solving 1D beam bending problems.
% This program is adapted from one provided on the Student Companion
% Website for "A First Course in Finite Elements" by Fish and Belytschko. 
%
% To change the problem being analyzed, edit the scripts BeamMesh and 
% InputData to create the proper geometry and boundary conditions. 
% 
% All variables are stored in structures, or "structs" which are detailed
% here:
% meshStruct.nCoords   % global nodal coordinates
% meshStruct.elCon     % connectivity array
% meshStruct.nnpe      % number of nodes per element (2)
% meshStruct.numDim    % number of spatial dimensions (1)
% meshStruct.numDOF    % degrees of freedom per node
% meshStruct.numNodes  % number of nodes
% meshStruct.numEq     % number of equations in global system
% meshStruct.numEls    % number of elements
% meshStruct.numSpans  % number of spans between points of interest in the
                       % beam
% meshStruct.gatherMat % global DOF for local DOF
% meshStruct.spanNum   % span number for each element
% meshStruct.elEI      % cross-section property for each element
% meshStruct.elDistLoad% distributed load for each element
% 
% globalSystem.K           % global stiffness matrix
% globalSystem.F           % global force vector including reaction vector
% globalSystem.d           % global solution vector including essential BCs
% globalSystem.reactionVec % reaction vector on essential DOF
% globalSystem.strain      % strain in each element
% globalSystem.stress      % stress in each element
% globalSystem.force       % internal force in each element
% 
% boundStruct.essBCs     % information about essential BCs
% boundStruct.appForces  % information about applied forces


% last edit: 27 September 2017 H. Ritz

% Preliminary steps
home; clear; close all;  % clean the workspace

% Preprocessing
meshStruct=BeamMesh;       % make the mesh geometry
[globalSystem,boundStruct,meshStruct]=...
    InputData(meshStruct); % read problem data such as 
                           % beam properties and BCs 

% Calculation and assembly of global stiffness matrix
globalSystem=...
    Assembly(globalSystem,meshStruct);  % creates local stiffness matrices 
                                        % and uses connectivity array to
                                        % create global stiffness matrix. 
                                        % This function takes advantage of 
                                        % the MATLAB 'sparse' command

% Solution Phase
globalSystem=...
    Soln(globalSystem,meshStruct,boundStruct); % applies essential BCs and 
                                               % solves the global system 
                                               % for nodal displacements 
                                               % and reaction forces

% Postprocessor Phase 
globalSystem=...
    PostProcess(globalSystem,meshStruct,boundStruct); % calculate stresses, 
                                                      % show deformed 
                                                      % shape, etc.
