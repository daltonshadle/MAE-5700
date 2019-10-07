function meshStruct=BeamMesh;
% meshStruct=BEAMMESH 
% This file defines the geometry for BEAM finite element code. You
% must define the "points of interest" for your problem and the number of
% elements to use for each span. Points of interest include the ends of the
% beam, the start and end locations for any distributed loads, the
% locations of any essential BCs or applied loads, etc. The POI should be
% given in order from left to right. A span is the region between
% consecutive POI. 

% last edit: 5 August 2015 H. Ritz

% some information about the types of elements
nnpe=2;   % number of nodes per element. 
numDOF=2; % number of degrees of freedom per node.
numDim=1; % number of spatial dimensions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT THIS FOR EACH PROBLEM
% Define the points of interest on the beam. This should include the left
% end and the right end as well as the locations of any applied loads or
% BCs and the beginning and end of any distributed loads.
pointsOfInterest=[0 4 8 12];
% For each span (i.e. the lengths between consecutive points of interest)
% define how many elements you want the span to have. 
% spanEls(i)=number of elements between POI #i and POI #i+1
spanEls=[1 1 1];
% Also define an alternative maximum element length. The code will either
% use this value or the number of elements per span, whichever creates more
% elements.
maxLength=12;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% use the points of interest and span element information to find the nodal
% point locations and make the connectivity array.
numPts=length(pointsOfInterest); % count how many points of interest
numSpans=length(spanEls);        % count how many spans
if numPts~=numSpans+1 % check that the POI and spans agree with each other
    error('Wrong number of spans indicated.')
end
nCoords=[];%initialize nodal vector
spanNum=[];%vector to keep track of which span each element is in
for sp=1:numSpans % loop over the spans
    pt2=pointsOfInterest(sp+1); % find the endpoints
    pt1=pointsOfInterest(sp);   % find the endpoints
    L=pt2-pt1;                  % span length
    spanElLeng=L/spanEls(sp);   % length of each element in the span if 
                                % using requested number of elements 
    % check whether we need to put more elements to satisfy maximum element
    % length and define the number of elements in this span accordingly
    if spanElLeng<maxLength     
        numSpanEls=spanEls(sp);  
    else
        numSpanEls=ceil(L/maxLength);
    end
    
    % now evenly space nodes along this span
    spanNodes=linspace(pt1,pt2,numSpanEls+1); 
    % add these new nodes to the end of the already existing nodal vector
    nCoords=[nCoords, spanNodes]; 
    spanNum=[spanNum; sp*ones(numSpanEls,1)]; % ditto for span number
end
% some of the nodes have been generated twice. (these correspond to the
% POI.) remove the duplicates.
nCoords=unique(nCoords)'; % keep only unique nodes
numNodes=size(nCoords,1); % count the number of nodes
elCon=[1:(numNodes-1); 2:numNodes]'; % make the element connectivity array
numEls=size(elCon,1);     % count the number of elements

% Use the connectivity array to define gatherMat which shows the global
% degrees of freedom for each local degree of freedom.
gatherMat=zeros(numEls,(nnpe*numDOF));
for n=1:nnpe  % loop over the number of nodes per element
    globalNodes=elCon(:,n); % global node number for this local node
    for d=1:numDOF % loop over the number of degrees of freedom per node
        % use global node numbers to find global DOFs
        % corresponding to the local DOFs.
        gatherMat(:,(n-1)*numDOF+d)=(globalNodes-1)*numDOF+d;
    end
end

% Package variables into the mesh struct
meshStruct.nCoords  =nCoords;  % global nodal coordinates
meshStruct.elCon    =elCon;    % connectivity array
meshStruct.nnpe     =nnpe;     % number of nodes per element (2)
meshStruct.numDim   =numDim;   % number of spatial dimensions (1)
meshStruct.numDOF   =numDOF;   % degrees of freedom per node
meshStruct.numNodes =numNodes; % number of nodes
meshStruct.numEls   =numEls;   % number of elements
meshStruct.gatherMat=gatherMat;% global DOF for local DOF
meshStruct.spanNum  = spanNum; % span number for each element
meshStruct.numSpans = numSpans;% number of spans in the mesh
meshStruct.pointsOfInterest = pointsOfInterest; % span endpoints