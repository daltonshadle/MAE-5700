function [globalSystem,boundStruct,meshStruct]=InputData(meshStruct);
% [globalSystem,boundStruct,meshStruct]=INPUTDATA(meshStruct)
% This file defines the beam properties, essential boundary conditions, and
% applied loads (both point and distributed) for the BEAM code. 
% The portions of the code that might change for each new
% problem are clearly indicated. 
% last edit: 5 August 2015 H. Ritz 

% unpack necessary input
numSpans=meshStruct.numSpans;
spanNum =meshStruct.spanNum;
pointsOfInterest=meshStruct.pointsOfInterest;
elCon   =meshStruct.elCon;
nCoords =meshStruct.nCoords;
numDOF  =meshStruct.numDOF;
numNodes=meshStruct.numNodes;
 
% Span properties (be sure to use consistent units). These properties
% may be the same or different for each span in the mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE THIS FOR EACH PROBLEM
spanEI 	= 1e4*ones(numSpans,1);   	% beam cross-section property 
spanDistributedLoad = [-1 -1 0]';   % magnitude of constant distributed 
                                    % transverse load on each span
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% spread the span properties to the elements
elEI=spanEI(spanNum);
elDistLoad=spanDistributedLoad(spanNum);


% Prescribed displacement or slope boundary conditions. 
% Each row is the number of the point of interest (numbered from left to
% right), the DOF, and the value for any essential BCs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE THIS FOR EACH PROBLEM
essBCs=[1 1 0; 
        1 2 0]; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Applied forces and moments. Each row is the number of the point of
% interest (numbered from left to right), the DOF, and the value for any
% applied loads.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE THIS FOR EACH PROBLEM
appForces=[2 1 -10; % for example, appForces=[3 2 20;4 1 -10;] means that 
           4 1 -20; % at the third point of interest there is an 
           4 2 20;  % applied moment with magnitude 20 in the CCW direction
           3 1 5];  % and at the fourth POI there is an applied transverse 
                    % load downward with magnitude 10.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize global system of equations
numEq=numNodes*numDOF;
F=zeros(numEq,1);
d=zeros(numEq,1);
K=zeros(numEq);

% Map the applied loads to the proper location in the global force vector
for frc=1:size(appForces,1)
    poi=pointsOfInterest(appForces(frc,1));  % POI for this applied load
    gnn=find(nCoords==poi);% global node number for this applied load
    gdof=(gnn-1)*numDOF+appForces(frc,2); % global DOF for the applied load
    val=appForces(frc,3);  % value of the applied load

    F(gdof)=val; % populate the global force vector with applied loads
end

% Package variables into the output structs
globalSystem.K=K;
globalSystem.F=F;
globalSystem.d=d;

boundStruct.essBCs    =essBCs;
boundStruct.appForces =appForces;

meshStruct.numEq =numEq;
meshStruct.elEI  =elEI;
meshStruct.elDistLoad=elDistLoad;

 