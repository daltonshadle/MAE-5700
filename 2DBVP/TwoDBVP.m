%% TWODBVP
%
% This code solves 2D boundary value problems of the form
% -grad . ([k]gradu) =f(x,y) using 3-node
% or 4-node elements. It is heavily based on N. Zabaras's code.
%
% The list of functions which must be edited for each problem is:
% MeshGeometry      : define the domain, element type, set plotting 
%                     options, plot the mesh, etc.
% InputData         : initialize BC types, define natural and essential BC
%                     values.
% DD                : constitutive relation
% FF                : define f(x,y) [RHS of BVP]
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
%
% boundStruct : defined in BoxGrid_2D (or loadFromGridFile) and InputData
%     boundStruct.elements
%         boundStruct.elements.Elems            : element numbers
%                                                 on each boundary
%         boundStruct.elements.SurfaceIndicator : flag for which element
%                                                 edge is on boundary
%     boundStruct.nodes
%         boundStruct.nodes.Nodes               : global node number on
%                                                 each boundary
%     boundStruct.SurfNat                       : prescribed natural BC
%     boundStruct.SurfEssV                      : prescribed essential BC
% 
% PlotInstructions : defined in MeshGeometry
%     PlotInstructions.plot_mesh     : flag for plotting mesh
%     PlotInstructions.plot_node     : flag for plotting node numbers
%     PlotInstructions.plot_boundary : flag for plotting boundary numbers
%     PlotInstructions.plot_contour  : flag for plotting nodal solution
%     PlotInstructions.plot_vector   : flag for plotting gradient (flux)
%
% globalSystem: defined in Assembly
%     globalSystem.K           % global stiffness matrix
%     globalSystem.F           % global force vector including reaction vector
%     globalSystem.d           % global solution vector including essential BCs
%     globalSystem.reactionVec % reaction vector on essential DOF

% last update: 7 November 2017 H. Ritz

% preliminary necessities
home; clear; close all; 
tic % start a timer, just for fun

% set up the mesh and define the boundaries
[meshStruct,boundStruct,PlotInstructions] = MeshGeometry; 
% choose plotting options and define BC type
boundStruct = InputData(boundStruct);
toc
disp(sprintf('\b   (Preprocessing)')) % output the time for meshing
tic
globalSystem = Assembly(meshStruct);% Loop over all the elements and 
                           % assemble the global sparse format of the 
                           % "stiffness" matrix and "force" vector                
toc
disp(sprintf('\b   (Assembly)')) % output the time for assembly
tic
% Apply different boundary conditions 
[globalSystem,boundStruct] = ApplyBC(boundStruct,meshStruct,globalSystem); 

% Solve the global system
globalSystem = Soln(globalSystem,boundStruct); 
toc
disp(sprintf('\b   (Solution)')) % output the time for the solution
tic
% Post-process for plots, flux, etc.
PostProcessor(PlotInstructions,meshStruct,globalSystem);
toc
disp(sprintf('\b   (Postprocessing)')) % output time for post-processing
