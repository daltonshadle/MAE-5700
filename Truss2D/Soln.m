function globalSystem = Soln(globalSystem,meshStruct,boundStruct)
% globalSystem = SOLN(globalSystem,meshStruct,boundStruct)
% Apply the essential boundary conditions and solve the global system for
% the nodal displacements for TRUSS2D. 
% This is solving the system [K]{d}={F}-{R}
% last edit: 10 July 2015 H. Ritz

% unpack necessary input
K=globalSystem.K;
F=globalSystem.F;
d=globalSystem.d;
essBCs=boundStruct.essBCs;
numDOF=meshStruct.numDOF;
numEq =meshStruct.numEq;

% unpack multipoint constraint variables
c = globalSystem.c;
q = globalSystem.q;
epsilon = globalSystem.epsilon; 
MPC_DOF = globalSystem.globMPC;

% assign variables and operations for multipoint constraints
% epsilonCTC is the matrix define by epsilon * C^T * C for calculations
epsilonCTC = epsilon .* c' * c;

% Add epsilonCTC to necessary parts of the equation
K = K + epsilonCTC;
globalSystem.K = K;

% partition the matrix K, vectors f and d
% first we want to find the global DOFs with essential boundary conditions
numEBC=size(essBCs,1); % number of essential boundary conditions
essDOF=zeros(numEBC,1);% initialize array of global indices to essential 
                       % boundary conditions
for ebc=1:numEBC
    % essDOF stores the index to the degrees of freedom with 
    % essential boundary conditions
   essDOF(ebc)=(essBCs(ebc,1)-1)*numDOF+essBCs(ebc,2);  
end
indF=setdiff(1:numEq,essDOF); % this returns the indices to the DOF that 
                              % do NOT have essential boundary conditions 
                              % (free DOF)
K_E	= K(essDOF,essDOF);       % Extract K_E matrix 
K_F	= K(indF,indF);           % Extract K_F matrix 
K_EF    = K(essDOF,indF);     % Extract K_EF matrix
f_F  	= F(indF);            % Extract f_F vector
d_E  	= essBCs(:,3);        % Extract d_E vector
 
forceMPC = epsilon*c'*q;
f_q = forceMPC(indF);
f_qe = forceMPC(essDOF);
% solve for d_F
d_F	=K_F\( f_F + f_q  - K_EF'* d_E);
 
% reconstruct the global solution vector
d(essDOF)=d_E;                
d(indF)=d_F;

% compute the reaction forces on the DOF with essential BCs
reactionVec = K_E*d_E + K_EF*d_F - F(essDOF) - f_qe;
% Package variables into the output structs
globalSystem.d=d;
globalSystem.reactionVec=reactionVec;



