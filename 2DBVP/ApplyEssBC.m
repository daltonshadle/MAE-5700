% ------------------------------------------------------------------------|
%                                                                         |
% MAE4700-5700, Finite Element Analysis for Mechanical & Aerospace Design |
%                                                                         |
% Copyright: Cornell University (this software should not be used without |
% written permission)                                                     |
%                                                                         |
% Authors: N. Zabaras (zabaras@cornell.edu) & Xiang Ma (xm25@cornell.edu) |
% Last update: 9 November 2015 Y. Xu                                      |                                  |
% ------------------------------------------------------------------------|
%
function  [bc, vals]=ApplyEssBC(i,boundStruct)

% It fills essential boundary conditions on boundary sideInd.
%
% Input,  sideInd      -- The boundary indicator/marker
%
% Output, bc           -- The list of degree of freedom of the essential
%                         boundary conditions at this node (just the node
%                         value for scalar fields)
%
% Output, vals         -- The list of values of essential boundary
%                         conditions at this node (here scalar)

%  Extract the node number of the boundary nodes on boundary ID.
bn = boundStruct.SurfEssV(i,1);

bc = boundStruct.nodes(bn).Nodes; % This is only valid for scaler field problem
Num = length(bc);
vals = ones(Num,1)*boundStruct.SurfEssV(i,2);




