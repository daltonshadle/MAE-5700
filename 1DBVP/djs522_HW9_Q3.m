% *************************************************************************
% MAE-5700 HW9 Q3
% Authors: Dalton Shadle
% Date: 11/6/19
% Notes: Lines of code listed here are the changes made to the original
%        1DBVP package, only the changes are listed.
% *************************************************************************

% From BarMesh.m ********************************************************
% ...
case 3
    % for quad elements
    numNodes=2*numEls+1;
    
    % coordinates of the evenly spaced nodes
    nCoords=linspace(xMin,xMax,numNodes);
    
    % element connectivity array, need connection of every node in element
    elCon=[1:2:numNodes-2; 2:2:numNodes-1; 3:2:numNodes]'; 
case 4
    % for cubic elements
    numNodes=3*numEls+1; 
    
    % coordinates of the evenly spaced nodes
    nCoords=linspace(xMin,xMax,numNodes); 
    
    % element connectivity array, need connection of every node in element
    elCon=[1:3:numNodes-3; 2:3:numNodes-2; 3:3:numNodes-1; 4:3:numNodes]';
    
    
% From Nmatrix.m ********************************************************
% ...
case 3 % quadratic elements
    N=[(1/2)*xi*(xi-1), -(1+xi)*(-1+xi), (1/2)*xi*(xi+1)];
case 4 % cubic elements
    N = [(9/16)*(1-xi)*(xi^2- 1/9), (27/16)*(xi^2-1)*(xi-1/3), ...
        (27/16)*(1-xi^2)*(xi+1/3), (9/16)*(1+xi)*(xi^2- 1/9)];

    
% From Bmatrix.m ********************************************************
% ...
case 3 % quadratic elements
    B=[xi-(1/2), -2*xi, xi+(1/2)];
case 4 % cubic elements
    B = [(9/16)*(-3*xi^2 + 2*xi + 1/9), (27/16)*(3*xi^2 - 2*xi/3 - 1), ...
        (27/16)*(-3*xi^2 - 2*xi/3 + 1), (9/16)*(3*xi^2 + 2*xi - 1/9)];








