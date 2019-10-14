function [globalSystem,boundStruct,meshStruct]=InputData(meshStruct);
% [globalSystem,boundStruct,meshStruct]=INPUTDATA(meshStruct)
% This file defines the element properties, essential boundary conditions, 
% and applied loads (both point and distributed) for the FRAME code. 
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
spanEI 	= 2.5e7*ones(numSpans,1);   % beam cross-section property 
spanEA  = 2.5e9*ones(numSpans,1);   % beam cross-section property 
spanDistributedLoad = [-10e3 0 0]'; % magnitude of constant distributed 
                                    % transverse load on each span
                                    % note! the positive direction is
                                    % determined by the span connectivity
                                    % array from FrameMesh. span endpoint 1
                                    % is at x'=0 and span endpoint 2 is at
                                    % x'=+L. positive q points in the
                                    % positive y' direction.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% spread the span properties to the elements
elEI=spanEI(spanNum);
elEA=spanEA(spanNum);
elDistLoad=spanDistributedLoad(spanNum);

% Prescribed displacement boundary conditions. Each row is the number of
% the point of interest, the DOF, and the
% value for any essential BCs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE THIS FOR EACH PROBLEM
essBCs=[1 1 0 ;   % for example, essBCs=[3 2 0;] means that 
        1 2 0 ;   % POI number 3 has a required y displacement of 0
        3 2 0];   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Applied forces and moments. Each row is the number of the point of
% interest, the DOF, and the value for any applied loads.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE THIS FOR EACH PROBLEM
appForces=[]; % for example, appForces=[3 3 20;4 1 -10;] means that 
%            4 1 -20; % at the third point of interest there is an 
%            4 2 20;  % applied moment with magnitude 20 in the CCW direction
%            3 1 5];  % and at the fourth POI there is an applied -x direction 
                    % load with magnitude 10.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize global system of equations
numEq=numNodes*numDOF;
F=zeros(numEq,1);
d=zeros(numEq,1);
K=zeros(numEq);

% Map the applied loads to the proper location in the global force vector
for frc=1:size(appForces,1)
    poi=pointsOfInterest(appForces(frc,1),:);  % POI for this applied load
    error(['Find the global node number corresponding to this POI. Possibly helpful commands include FIND and ISMEMBER.'])
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
meshStruct.elEA  =elEA;
meshStruct.elDistLoad=elDistLoad;



 