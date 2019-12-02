function  [essDOF, essVals]=ApplyEssBC(i,boundStruct)
% ApplyEssBC deals with essential boundary conditions and outputs a column
% vector essDOF--the list of degree of freedom of the essential BCs along
% the current surface, and a column vector essVals--the corresponding 
% values of the essential BCs.

% last update: 17 Nov 2015 H. Ritz; Y. Xu 

%  Extract the node number of the boundary nodes on boundary ID.
bn = boundStruct.SurfEssV(i,1); % current surface number
dir = boundStruct.SurfEssV(i,2); % direction (x or y)
BNodes = boundStruct.nodes(bn).Nodes; % Boundary node numbers

essDOF = 2*(BNodes-1) + dir; % degree of freedom
essDOF = essDOF';

essVals = ones(length(essDOF),1)*boundStruct.SurfEssV(i,3); % or use repmat function 