function [meshStruct,boundStruct]=InputData(meshStruct,boundStruct)
% [meshStruct,boundStruct]=InputData(meshStruct,boundStruct);
% Define the essential and natural BCs, and define the mateiral stiffness
% matrix D

% last update: 15 Nov 2015 H. Ritz; Y. Xu

% For complex geometries, you need to firstly plot the mesh (set both 
% PlotInstructions.plot_mesh and PlotInstructions.plot_boundary to be
% 'yes') to see how the number of boudaries are defined.

% Define the essential BCs
% boundStruct.SurfEssV = [];
boundStruct.SurfEssV = [5 1 0
                        4 2 0]; % e.g. [4 2 20] means all nodes on surface # 4,
%                                % degree of freedom #2 (y direction), has a value of 20.

% Define the natural BCs
% The natural boundary condition is defined in tangential and normal
% direction (rather than global x and y direction),outer normal is
% positive.
boundStruct.SurfNat = [3 0 -1e7]; % e.g. [3 10 -10] means surface # 3 has 
                                 % a constantly distributed tangential traction
                                 % 10 and normal traction (pointing in) 10.
% Define material properties
E          =2e11; % Young's Modulus
nu         =0.35; % Poisson's Ratio
PlaneStress='yes';% 'yes' for plane stress, 'no' for plane strain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch PlaneStress
    case 'yes'
        D=E/(1-nu^2)*[1 nu 0 ; nu 1 0 ; 0 0 (1-nu)/2];
    case 'no'
        D=E/((1+nu)*(1-2*nu))*[(1-nu) nu 0 ; nu (1-nu) 0 ; 0 0 (1-2*nu)/2];
    otherwise
        error('Is this plane stress or plane strain?');
end

meshStruct.Material.D=D;
meshStruct.Material.E=E;
meshStruct.Material.nu=nu;

