% *************************************************************************
% Filename: Uy_HW12_P2.m
% Authors: Dalton and Sairam
% Date: 11/26/19
% Notes: This script is for finding the Uy displacement of the beam at x=L
%     for HW12 Problem 2.
% *************************************************************************

% unpack variables from structures
d = globalSystem.d;
nCoords = meshStruct.nCoords;

% get y-displacements
y_dof_index = 2*[1:length(d)/2];
Uy_d = d(y_dof_index);

% find all nodes at x=L
L = max(nCoords(:,1));
L_nodes = find(nCoords(:,1) == L);

% get y-displacements for x=L nodes
Uy_L_nodes = Uy_d(L_nodes);
Uy_L = mean(Uy_L_nodes);

% print results
fprintf('U_y displacement at x=L: %0.6f\n', Uy_L);

fprintf('Max d: %0.6f\n', min(Uy_d));

% exact solution
E = 200e+9; % Pa
P_b = 50e+3; % kN/m
a = 0.1; % m
L = 1; % m

% Note: I = ba^3/12
exact_sol = - (P_b * L^3 * 12) / (3 * E * a^3);
fprintf('Exact Solution: %0.6f\n', exact_sol);

% Error calculations
L2_error = sqrt((Uy_L - exact_sol)^2) / sqrt(exact_sol^2);
fprintf('Error: %0.6f\n', L2_error);
