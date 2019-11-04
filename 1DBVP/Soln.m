function globalSystem = Soln(globalSystem,meshStruct,boundStruct)
% globalSystem = SOLN(globalSystem,meshStruct,boundStruct)
% Apply the essential boundary conditions and solve the global system for
% the nodal degrees of freedom for ONEDBVP. 
% last edit: 6 August 2015 H. Ritz

% unpack necessary input
K      =globalSystem.K;
F      =globalSystem.F;
essDOF =boundStruct.essDOF;
ebcVals=boundStruct.ebcVals;
numEq  =meshStruct.numEq;

% partition the matrix K, vectors f and d
freeDOF=setdiff(1:numEq,essDOF); % this returns the indices to the DOF that 
                             % do NOT have essential boundary conditions 
                             % (free DOF)

K_E	= K(essDOF,essDOF);      % Extract K_E matrix 
K_F	= K(freeDOF,freeDOF);    % Extract K_F matrix 
K_EF= K(essDOF,freeDOF);     % Extract K_EF matrix
f_F = F(freeDOF);            % Extract f_F vector
d_E = ebcVals;               % Extract d_E vector
 
% solve for d_F
d_F	=K_F\( f_F - K_EF'* d_E);
 
% reconstruct the global solution d
d=zeros(numEq,1);
d(essDOF)=d_E;                
d(freeDOF)=d_F;

% compute the reactions on the DOF with essential BCs
reactionVec = K_E*d_E+K_EF*d_F-F(essDOF);

% Package variables into the output struct
globalSystem.d=d;
globalSystem.reactionVec=reactionVec;

