function [meshStruct,boundStruct]=InputData(meshStruct,boundStruct)
% [meshStruct,boundStruct]=InputData(meshStruct,boundStruct);
% Define the essential and natural BCs, and define the mateiral stiffness
% matrix D

% last update: 15 Nov 2015 H. Ritz; Y. Xu

% For complex geometries, you need to firstly plot the mesh (set both 
% PlotInstructions.plot_mesh and PlotInstructions.plot_boundary to be
% 'yes') to see how the number of boudaries are defined.

%% Define the essential BCs
boundStruct.SurfEssV = [1 1 0; 1 2 0; 1 3 0;
                        2 1 0; 2 2 0; 2 3 0;
                        3 1 0; 3 2 0; 3 3 0;
                        4 1 0; 4 2 0; 4 3 0]; 
                    % [a b c] means on "a" surface, the "b" degree of
                    % freedom has a value of "c"
                    
%% Define the natural BCs
% The natural boundary condition is defined in tangential and normal
% direction (rather than global x and y direction),outer normal is
% positive.
boundStruct.SurfNat = []; % "Have to change"
%% Define material properties
E          =2e11; % Young's Modulus
nu         =0.35; % Poisson's Ratio
thickness = meshStruct.thickness; % uniform thickness over the plate
PlaneStress='yes';% 'yes' for plane stress, 'no' for plane strain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch PlaneStress
    case 'yes'
        D=E*thickness^3/(12*(1-nu^2))*[1 nu 0 ; nu 1 0 ; 0 0 (1-nu)/2];
    case 'no'
        D=E*thickness^3/(12*(1+nu)*(1-2*nu))*[(1-nu) nu 0 ; nu (1-nu) 0 ; 0 0 (1-2*nu)/2];
    otherwise
        error('Is this plane stress or plane strain?');
end
%% Package material properties

meshStruct.Material.D=D;
meshStruct.Material.E=E;
meshStruct.Material.nu=nu;

