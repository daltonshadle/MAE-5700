function [globalSystem,boundStruct,meshStruct]=InputData(meshStruct) 
% [globalSystem,boundStruct,meshStruct]=INPUTDATA(meshStruct)
% This file defines the element properties, prescribed
% displacements (essential boundary conditions) and applied loads for the
% TRUSS2D code. 
% The portions of the code that might change for each new
% problem are clearly indicated. 
% last edit: 10 July 2015 H. Ritz

 
% unpack necessary input
numEls = meshStruct.numEls;
numDOF = meshStruct.numDOF;
numNodes = meshStruct.numNodes;

% Element properties (be sure to use consistent units). These properties
% may be the same or different for each element in the mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE THIS FOR EACH PROBLEM
elArea 	= 1*10^(-4)*ones(numEls,1);   % Elements area, column vector where 
                                      % each component is the area of that 
                                      % element 
elYM    = 200*10^(9)*ones(numEls,1);  % Young's Modulus, column vector where 
                                      % each componenet is the YM of that 
                                      % element 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Applied forces. Each row is the global node number, the DOF, and
% the value for any applied loads
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE THIS FOR EACH PROBLEM
% better example [globalNodeNum, direction, force]
% where direction (x=1, y=2, z=3) and force unit 
% = Newtons
% for example, appForces=[3 2 20e3]; means that 
% global node number 3 has an applied load 
% in the y direction with magnitude 20e3
myForce=1000;
appForces=[2 1 myForce]; 

% For HW4
appForces=[3 2 -myForce]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prescribed displacement boundary conditions. Each row is the global node
% number, the DOF, and the value for any essential BCs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE THIS FOR EACH PROBLEM
% better example [globalNodeNum, direction, displacement]
% where direction (x=1, y=2, z=3) and displacement unit 
% = meters
% for example, essBCs=[3 2 0;] means that 
% global node number 3 has a required displacement  
% of 0 in the y direction
essBCs = [1 1 0;
          1 2 0];
      
% For HW4
essBCs = [1 1 0;
          1 2 0];

% variables for defining multipoint constraints
% epsilon (float) - penalty parameter for mp constraints
% multiConstraint (m x 2) - definition of multipoitn constraints where the
%     usage is on constraint per row as follows [node_#, theta_angle(radians)]
%     where node_# is the number of the node and theta_angle is the angle
%     of the sloped mp constraint in radians
epsilon = 1e15;
multiConstraint = [4, deg2rad(60)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize global system of equations
numEq = numNodes * numDOF;
F = zeros(numEq, 1);
d = zeros(numEq, 1);
K = zeros(numEq);
% initialize multipoint constraint variable and matrices
numMultiConstraint = size(multiConstraint, 1);
q = zeros(numMultiConstraint, 1);
c = zeros(numMultiConstraint, numEq);

% Map the applied loads to the proper location in the global force vector
for frc = 1:size(appForces,1)
    gnn = appForces(frc,1);        % global node number for this applied load
    gdof = (gnn-1) * numDOF + appForces(frc,2); % global DOF for the applied load
    val = appForces(frc,3);        % value of the applied load

    F(gdof) = val; % populate the global force vector with applied loads
end

% Package variables into the output structs
globalSystem.K = K;
globalSystem.F = F;
globalSystem.d = d;

% Add multipoint constraint variables to globalSystem struct
globalSystem.epsilon = epsilon;
globalSystem.c = c;
globalSystem.q = q;
globalSystem.multiConstraint = multiConstraint;

boundStruct.essBCs    = essBCs;
boundStruct.appForces = appForces;

meshStruct.numEq = numEq;
meshStruct.elArea= elArea;
meshStruct.elYM  = elYM;


 