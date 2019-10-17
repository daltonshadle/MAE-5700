function meshStruct=FrameMesh
% meshStruct=FRAMEMESH 
% This file defines the geometry for FRAME finite element code. You must
% define the "points of interest" for your problem and the number of
% elements to use for each span. Points of interest include the ends of the
% beam, the start and end locations for any distributed loads, the
% locations of any essential BCs or applied loads, etc. 
% A span connects POI.

% last edit: 5 August 2015 H. Ritz

% some information about the types of elements
nnpe=2;   % number of nodes per element. 
numDOF=3; % number of degrees of freedom per node.
numDim=2; % number of spatial dimensions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the points of interest on the frame. This should include the
% enpoints of all frame members as well as the locations of any applied
% loads or BCs and the beginning and end of any distributed loads. 
% First column is x coordinates, second column is y coordinates.
pointsOfInterest=[0 0;
                  0 3; 
                  0 6;
                  4 6;
                  4 0;
                  4 3];
pointsOfInterest=[0, 0;
                  0.5, sqrt(3)/2;
                  1 0];
% Define the spans and how many elements you want each span to have. A span
% can't have any POI inside it.
spanCon=[1 2;    % this is the connectivity array for the spans.
         2 3;
         2 6;
         3 4;
         4 6;
         6 5];
spanCon=[1, 2;
         2, 3;
         1, 3];
spanEls=[20 20 20 20 20 20]; % spanEls(S) is the number of elements span S should have
spanEls=[1, 1, 1];
% Also define an alternative maximum element length. The code will either
% use this value or the number of elements per span, whichever creates more
% elements.
maxLength=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% use the points of interest and span element information to find the nodal
% point locations and make the connectivity array.
numPts=length(pointsOfInterest); % count how many points of interest
numSpans=length(spanEls);        % count how many spans
nCoords=[];                      % initialize nodal coordinate array
elCon=[];                        % initialize element connectivity array
spanNum=[];                      % initialize vector to keep track of 
                                 % which span each element is in
for sp=1:numSpans % loop over the spans
    poi1=pointsOfInterest(spanCon(sp,1),:);% coordinates of 1st POI
    poi2=pointsOfInterest(spanCon(sp,2),:);% coordinates of 2nd POI
    x1=poi1(1); y1=poi1(2); % find the span endpoints
    x2=poi2(1); y2=poi2(2); % find the span endpoints
    L=sqrt((x2-x1)^2+(y2-y1)^2);% L = length of the span
    spanElLeng=L/spanEls(sp);   % length of each element in the span if 
                                % using requested number of elements 
    % check whether we need to put more elements to satisfy maximum element
    % length and define the number of elements in this span accordingly
    if spanElLeng<maxLength     
        numSpanEls=spanEls(sp); 
    else
        numSpanEls=ceil(L/maxLength);%
    end
    % find the x and y coordinates of evenly spaced nodes along this span
    spanX=linspace(x1,x2,numSpanEls+1); % x
    spanY=linspace(y1,y2,numSpanEls+1); % y
    nCoords=[nCoords; [spanX', spanY']];% add these new nodes to the end of
                                        % the already existing nodal array
    numNodes=size(nCoords,1);           % count how many nodes so far
    % make the element connectivity. This will have to be corrected when
    % we're deleting duplicate nodes, but for now we can connect the
    % elements in a line along the new nodes.
    elCon=[elCon;[(numNodes-numSpanEls):(numNodes-1); (numNodes-numSpanEls+1):(numNodes)]'];
    % keep track of the span number for each new element
    spanNum=[spanNum; sp*ones(numSpanEls,1)]; 
end
% there are several duplicate nodes. now we get rid of them and fix the
% connectivity array
% [nCoords,IA,IC]=unique(nCoords,'rows'); % choose the unique nodes 
[nCoords,IA,IC]=unique(nCoords,'rows','stable'); % choose the unique nodes 
                                                 % and save the indexes 
                                                 % to the proper rows
nCoords
elCon=IC(elCon);% fix the element connectivity array using 
                 % only the unique nodes
elCon              
numNodes=size(nCoords,1); % now save the final number of nodes 
numEls=size(elCon,1);     % and the final number of elements

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