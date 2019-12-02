% LinElast
%
% This code solves 2D linear elasticity problems using 3-node
% or 4-node elements. It is heavily based on N. Zabaras's code.
%
% The list of functions which must be edited for each problem is:
% MeshGeometry      : define the domain, element type, set plotting 
%                     options, plot the mesh, etc.
% InputData         : define natural and essential BC values, define
%                     material stiffness matrix D.
% BodyForce         : define x and y components of body force.
% 
% All of the structs used, including their fields and the
% functions in which they are first defined, are listed here:
% 
% meshStruct : defined in MeshGeometry
%     meshStruct.elCon     : connectivity array
%     meshStruct.gatherMat : gather/scatter matrix
%     meshStruct.nCoords   : nodal coordinates
%     meshStruct.nnpe      : number of nodes per element
%     meshStruct.nsd       : number of spatial dimensions
%     meshStruct.numDOF    : number of degrees of freedom per node
%     meshStruct.numEls    : number of elements
%     meshStruct.numEq     : number of equations in global system
%     meshStruct.numNodes  : number of global nodes
%     meshStruct.Material  : struct with material properties defined in
%                            InputData
%         meshStruct.Material.E : Young's modulus  
%         meshStruct.Material.nu: Poisson's ratio  
%         meshStruct.Material.D : constitutive relation matrix 
%                                 such that stress=D*strain  
%
% boundStruct : defined in BoxGrid_2D (or loadFromGridFile) and InputData
%     boundStruct.elements
%         boundStruct.elements.Elems            : element numbers
%                                                    on each boundary
%         boundStruct.elements.SurfaceIndicator : flag for which element
%                                                    edge is on boundary
%     boundStruct.nodes
%         boundStruct.nodes.Nodes               : global node number on
%                                                    each boundary
%     boundStruct.SurfNat                       : prescribed natural BC
%     boundStruct.SurfEssV                      : prescribed essential BC
% 
% PlotInstructions : defined in InputData
%     PlotInstructions.plot_mesh     : flag for plotting mesh
%     PlotInstructions.plot_node     : flag for plotting node numbers
%     PlotInstructions.plot_boundary : flag for plotting boundary numbers
%     PlotInstructions.plot_contour  : flag for plotting nodal solution
%     PlotInstructions.plot_deformed : flag for plotting deformed mesh
%     PlotInstructions.plot_fringes  : flag for plotting fringe pattern

% last update: 15 Nov 2015 H. Ritz; Y. Xu

%%  Clearing windows and variables

clc; clear variables; close all; 

%% set up the mesh and define the boundaries

[meshStruct,boundStruct,PlotInstructions]=MeshGeometry;

%%

%% Preprocessing
tic
[meshStruct,boundStruct]=InputData(meshStruct,boundStruct);
toc
fprintf('\b   (Preprocessing)\n') % output the time for meshing

%% Assembly
tic
globalSystem = Assembly(meshStruct);% Loop over all the elements and 
                                  % assemble the global sparse format of the 
                                  % "stiffness" matrix and "force" vector                
toc
fprintf('\b   (Assembly)\n') % output the time for assembly

%% Boundary conditions
tic
% Apply different boundary conditions 
[globalSystem,boundStruct] = ApplyBC(boundStruct,meshStruct,globalSystem); 
toc
fprintf('\b   (ApplyBC)\n') % output the time for assembly

%% Solution: Solve the global system
tic
globalSystem = Soln(globalSystem,boundStruct);
toc
fprintf('\b   (Solution)\n') % output the time for the solution

%% Post-process for strains, stresses and moments
tic
PostProcessor(PlotInstructions,meshStruct,globalSystem);
toc
fprintf('\b   (Postprocessing)\n') % output time for post-processing
%%
