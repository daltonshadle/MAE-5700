% *************************************************************************
% Name: InputData.m
% Notes: [meshStruct,boundStruct]=InputData(meshStruct,boundStruct);
%     Define the essential and natural BCs, and define the mateiral stiffness
%     matrix D
%     last update: 15 Nov 2015 H. Ritz; Y. Xu
% Project Updates: Pulled from LinElast code base from MAE-5700 course,
%     updates were made to implement Kirchhoff Theory Plate Bending
% Update Authors: Dalton and Sairam
% *************************************************************************

function [meshStruct,boundStruct]=InputData(meshStruct,boundStruct)
%% Define the essential BCs
% takes the form [a, b, c]
% - a = surface number
% - b = DOF (for plate bending 1=u_z, 2=theta_x, 3=theta_y)
% - c = value of BC

% Simply Supported Essential BC's
boundStruct.SurfEssV = [1 1 0; 2 1 0; 3 1 0; 4 1 0];

% Clamped Essential BC's
boundStruct.SurfEssV = [1 1 0; 1 2 0; 1 3 0;
                        2 1 0; 2 2 0; 2 3 0;
                        3 1 0; 3 2 0; 3 3 0;
                        4 1 0; 4 2 0; 4 3 0;];
                   
                    
%% Define the natural BCs
% takes the form [a, b, c]
% - a = surface number
% - b = force DOF (for plate bending 1=f_z, 2=M_x, 3=M_y)
% - c = value of BC
% NOTE: M_x refers to the moment applied along a particular x-plane
%       M_y refers to the moment applied along a particular y-plane
% SURFACE NOTE: Surface 1: Right-Hand Convention
%               Surface 2: Left-Hand Convention
%               Surface 3: Right-Hand Convention
%               Surface 4: Left-Hand Convention

% Simiply Supported Natural BC's
boundStruct.SurfNat = [1 2 0; 1 3 0;
                       2 2 0; 2 3 0;
                       3 2 0; 3 3 0;
                       4 2 0; 4 3 0];
                   
% Clamped Nautral BC's
boundStruct.SurfNat = [];

%% Define material properties
E          =2e11; % Young's Modulus
nu         =0.3; % Poisson's Ratio
thickness  = meshStruct.thickness; % uniform thickness over the plate
PlaneStress='yes';% 'yes' for plane stress, 'no' for plane strain

% Note: for Kirchhoff Theory, this will always be plane stress
switch PlaneStress
    case 'yes'
        D=E/((1-nu^2))*[1 nu 0 ; nu 1 0 ; 0 0 (1-nu)/2];
    case 'no'
        D=E/((1+nu)*(1-2*nu))*[(1-nu) nu 0 ; nu (1-nu) 0 ; 0 0 (1-2*nu)/2];
    otherwise
        error('Is this plane stress or plane strain?');
end

%% Package material properties
meshStruct.Material.D=D;
meshStruct.Material.E=E;
meshStruct.Material.nu=nu;

