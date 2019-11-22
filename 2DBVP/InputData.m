function [boundStruct]=InputData(boundStruct)
% [boundStruct]=InputData(boundStruct);
% Define the essential and natural BC.
% Set this information for each problem
% last update: 7 November 2017 H. Ritz

% For complex geometries, you need to first plot the mesh (set both 
% PlotInstructions.plot_mesh and PlotInstructions.plot_boundary to be
% 'yes') to see how the number of boudaries are defined.

% boundStruct.SurfEssV = [3 1000; ]; % e.g. [3 20] means all nodes on surface # 3 
%                                    % have a fixed value of 20.

% boundStruct.SurfEssV = [3 100 ; 1 0 ; 2 0]; % this is for the example of plate with a hole 
boundStruct.SurfEssV = [1 0]; % this is for torsion of quarter circle shaft

% boundStruct.SurfEssV = [1 0; 2 0; 3 0; 4 0]; % for rectangle

% The natural boundary condition is expressed as
% k*(gradu) . n = alpha * u + beta
% for Neumann BC, alpha = 0; beta = - q_bar
% alpha is nonzero for convection BC (mixed BC).

% BE CAREFUL ABOUT THE SIGN CONVENTIONS HERE
% boundStruct.SurfNat = [1 0 0; 2 0 0]; % e.g. [2 5 -10] means  surface # 2 has 
                                 % natrual BC defined as:
                                 %  k*(gradu) . n = 5 * u -10
                                 % k is defined in the function DD.m, and
                                 % can generally depend on positions.
% boundStruct.SurfNat = [4 0 0; 5 0 0]; % this is for the example of plate with a hole 
boundStruct.SurfNat = [2 0 0; 3 0 0]; % this is for torsion of quarter circle shaft
% boundStruct.SurfNat = []; % for rectangle


% point source boundary condition for heat convection
% in the for [heat source, x-coord, y-coord]
% boundStruct.PointSource = [250 0 5e-2];

