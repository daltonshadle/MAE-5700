function [boundStruct]=InputData(meshStruct);
% [boundStruct]=INPUTDATA(meshStruct);
% This file defines the boundry conditions and point forces for the 1D
% boundary value problem. Unlike the truss code, no material properties are
% set here. This is because this code can be used to solve several
% different problems, such as linear elasticity or heat conduction, or any
% problem whose strong form can be written: 
% -d/dx( p(x)du/dx )=f(x)
% Since the problem statement is general, any material properties belong in
% the scripts PP, or FF. The portions of the code that might change for
% each new problem are clearly indicated.
% last edit: 5 August 2015 H. Ritz 

% unpack necessary variables

% Applied "forces". Each row has the x position of any point "force," and
% the value of the point "force."
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE THIS FOR EACH PROBLEM
appForces=[];
% appForces=[5 24]; % for example, appForces=[0 10; 5 -24] means that 
%             % at x=0 there is an applied point force of 10 and 
%             % at x=5 there is an applied point force of -24.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Boundary conditions. In the 1D problem you must define 2 boundary
% conditions: one for the left side of the domain and one for the right
% side of the domain.
% boundCond=[T1 V1;
%            T2 V2];
% where T1 is the flag for type of boundary condition on the left side of
% the domain, V1 is the value of the boundary condition on the left side of
% the domain, and T2 and V2 are the type and value of the boundary
% conditions on the right side of the domain. The type T1 and T2 should be 
% 1 for an essential boundary condition (given u) and 
% 2 for a natural boundary condition (given du/dx).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE THIS FOR EACH PROBLEM
boundCond=[1 0;       
           2 4*pi*cos(4*pi)];   
% boundCond=[1 0;       
%            2 0];   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Package output structs
boundStruct.boundCond=boundCond;
boundStruct.appForces=appForces;

 