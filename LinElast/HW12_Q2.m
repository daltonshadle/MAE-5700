% *************************************************************************
% Name: HW12_Q2.m
% Authors: Dalton and Sairam
% Date: 12/7/19
% Notes: Summarized code changes for HW12 Question 2
% *************************************************************************


% From InputData.m *******************************************************
% ...
% Essential BC's for HW12 Problem 2, fixed on surface 4
boundStruct.SurfEssV = [4 1 0;
                        4 2 0];
% ...
yt = meshStruct.yt;
% Natural BC's for HW12 Problem 2, traction free on surface 1,3, transverse
% traction on surface 2 equal to P/(b*a)
boundStruct.SurfNat = [2 -50e+3/yt 0;
                       1 0 0;
                       3 0 0];

% From MeshGeomtry.m ******************************************************
% ...
xl = 0.0;           % left location of the range in the x direction
xr = 1.0;           % right location of the range in the x direction
yb = 0.0;           % bottom location of the range in the y direction
yt = 0.1;           % top location of the range in the y direction
% NOTE: yt = 0.01 for part (c)


% From Uy_HW12_P2.m *******************************************************
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
fprintf('U_y displacement at x=L: %0.8f\n', Uy_L*1000); % in mm

% exact solution
E = 200e+9; % Pa
P_b = 50e+3; % kN/m
a = 0.01; % m
L = 1; % m

% Note: I = ba^3/12
% Solving: u_y = P*L^3/(3*E*I)
exact_sol = - (P_b * L^3 * 12) / (3 * E * a^3);
fprintf('Exact Solution: %0.8f\n', exact_sol*1000); % in mm

% L2_Error calculations
% NOTE: since U_y is constant at x=L in the FEM formulation and is also
% constant in the exact solution, the error is calculated as such
L2_error = sqrt((Uy_L - exact_sol)^2) / sqrt(exact_sol^2);
fprintf('Error: %0.6f\n', L2_error);













