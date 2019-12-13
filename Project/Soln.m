% SOLN 
% Apply the essential boundary conditions and solve the global system for
% the nodal degrees of freedom for LinElast. 
% last edit: 22 November 2017 H. Ritz

function globalSystem = Soln(globalSystem,boundStruct)
% unpack the things you need
K = globalSystem.K;
F = globalSystem.F;
essDOF = boundStruct.essDOF;
ebcVals = boundStruct.ebcVals;

numEq=length(F);
% partition the matrix K, vectors f and d
freeDOF=setdiff(1:numEq,essDOF); % this returns the indices to the DOF that 
                             % do NOT have essential boundary conditions 
                             % (free DOF)

K_E	= K(essDOF,essDOF);      % Extract K_E matrix 
K_F	= K(freeDOF,freeDOF);    % Extract K_F matrix 
K_EF= K(essDOF,freeDOF);     % Extract K_EF matrix
f_F = F(freeDOF);           % Extract f_F vector
d_E = ebcVals;               % Extract d_E vector
 
% solve for d_F
if isempty(d_E)
    d_F	=K_F\( f_F );
else
    d_F	=K_F\( f_F - K_EF'* d_E);
end
 
% reconstruct the global solution d
d(essDOF)=d_E;                
d(freeDOF)=d_F;

% compute the reaction forces on the DOF with essential BCs
if isempty(d_E)
    reactionVec = K_EF*d_F;
else
    reactionVec = K_E*d_E+K_EF*d_F-F(essDOF);
end

globalSystem.d = d;
globalSystem.reactionVec = reactionVec;


