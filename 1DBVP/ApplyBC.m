function [globalSystem,boundStruct]=ApplyBC(globalSystem, meshStruct, boundStruct)
% globalSystem=APPLYBC(globalSystem, meshStruct, boundStruct)
% This function applies both the essential and natural boundary conditions.
% Essential boundary conditions are of the form u(x)=Value
% Natural boundary conditions are of the form du(x)/dx=Value
% 1D problems will have exactly two boundary conditions (one on the left
% one on the right)

% last edit: 5 August 2015 H. Ritz

% unpack necessary input
F        =globalSystem.F;
nCoords  =meshStruct.nCoords;
numNodes =meshStruct.numNodes;
boundCond=boundStruct.boundCond;

% unpack the boundCond array
BCtype=boundCond(:,1); % first column indicates the type of BC
BCval =boundCond(:,2); % second column stores the value of the BC
BCnode=[1 numNodes];

essDOF=[]; ebcVals=[];  % Initialize debc, ebcVals
% There are two boundary nodes for 1D problem
for bn = 1:2 % for each boundary node...
    switch BCtype(bn)
        case 1 % essential BC
            essDOF =  [essDOF BCnode(bn)];  % this essential BC information 
            ebcVals = [ebcVals; BCval(bn)]; % is used in Soln script
        case 2 % natural BC
            val=BCval(bn);
            if (bn==1) % left boundary needs the sign fixed
                val=-val;  
            end
            x = nCoords(BCnode(bn));
            % natural boundary condition changes global force vector
            F(BCnode(bn)) = F(BCnode(bn)) + PP(x)*val;          
        otherwise
            error('Invalid boundary condition code')
    end
end

% Package output structs
globalSystem.F=F;
boundStruct.essDOF=essDOF;
boundStruct.ebcVals=ebcVals;