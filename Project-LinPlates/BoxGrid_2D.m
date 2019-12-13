% ------------------------------------------------------------------------|
%                                                                         |
% MAE4700-5700, Finite Element Analysis for Mechanical & Aerospace Design |
%                                                                         |
% Copyright: Cornell University (this software should not be used without |
% written permission)                                                     |
%                                                                         |
% Authors: N. Zabaras (zabaras@cornell.edu) & Xiang Ma (xm25@cornell.edu) |
%                                                                         |
% ------------------------------------------------------------------------|
% You should not need to edit this function.
function [nodalcoords,connectivity,BoundaryStruct]=BoxGrid_2D(nsd, xl, xr, yb, yt, nnpe, nelemx, nelemy );

%  BoxGrid_2D produces a grid of pairs of 3 node triangles or 4-node
%  quadrilaterals in the rectangular domain [xl,xr] x [yb,yt].
%
%  For each element, the nodes and boundary indicators are listed in
%  counter-clockwise order. Note that for nelemx and nelemy are
%  defining the subdivisions in the x and y directions rather than
%  the number of elements -- see Figures below. 
%
%  Example:
%
%    Input:
%
%      NELEMX = 3, NELEMY = 2
%
%
%          Grid T3:			          Grid Q4:
%          side 3                      side 3 
%    9---10---11---12		      9---10---11---12
%    |\ 8 |\10 |\12 |		      |    |    |    |
%    | \  | \  | \  |		      |    |    |    |
%  s |  \ |  \ |  \ | s		    s |  4 |  5 |  6 | s
%  i |  7\|  9\| 11\| i		    i |    |    |    | i
%  d 5----6----7----8 d	        d 5----6----7----8 d
%  e |\ 2 |\ 4 |\ 6 | e	        e |    |    |    | e
%  4 | \  | \  | \  | 2		    4 |    |    |    | 2
%    |  \ |  \ |  \ |		      |  1 |  2 |  3 |
%    |  1\|  3\|  5\|		      |    |    |    |
%    1----2----3----4		      1----2----3----4
%         side 1                        side 1
%  Reference Element T3:	    Reference Element Q4:
%
%    |				               |           
%    1    3			               1    4--------3
%    |    |\			           |    |        |
%    |    | \			           |    |        |
%    Eta  |  \			           Eta  |        | 
%    |    |   \			           |    |        |
%    |    |    \		           |    |        |
%    -1   1-----2		          -1    1--------2
%    |				               |          
%    +-- -1--Xi--1-->		       +-- -1--Xi--1-->
%
%
%
%  Node labeling:
%
%    NW--NE
%     |\ |
%     | \|
%    SW--SE
%

%  Parameters:
%
%    Input, integer NELEMX, NELEMY, the number of subdivisions along the
%    X and Y directions.  The number of elements generated will be
%    NELEMX * NELEMY (for Q4 elements) and 2* NELEMX * NELEMY (for T3).
%
%    Input,  double xl, xr, the left and right range along X direction.
%    Input,  double yb, yt, the bottom and top range along Y direction.
%

switch nnpe
    case 4
        nel = nelemx*nelemy;                % number of total elements
    case 3
        nel = 2*nelemx*nelemy;
    otherwise
        error('This element type is not implemented\n');
end

nno = (nelemx+1)*(nelemy+1);            % number of total nodes
                              
nodalcoords = zeros(nno,nsd);           % initialize arrays
connectivity = zeros(nnpe,nel);

                                        % collect the coordinates of the nodes
 x = linspace(xl,xr,nelemx+1);          % equal distribution of the x nodes
 y = linspace(yb,yt,nelemy+1);          % equal distribution of the y nodes 

for j = 0 : nelemy                      % loop over all nodes, set up the coordinates 
    for i = 0 : nelemx

        nodeID = (nelemx + 1) * j  + i + 1 ;

        nodalcoords(nodeID,1) = x(i+1);       % Set up the coordinates
        nodalcoords(nodeID,2) = y(j+1);       
    end
end

                                        % What global nodes are in each
                                        % boundary segment?
                                  
BoundaryNodes(1).Nodes = [1:nelemx+1];  % The first boundary (bottom)

                                        % The second boundary (right)
BoundaryNodes(2).Nodes = [1:(nelemy+1)]*(nelemx+1);

                                        % The third boundary (top)
BoundaryNodes(3).Nodes = ( nelemx + 1) * nelemy  + [1 : nelemx+1];

                                        % The fourth boundary (left)
BoundaryNodes(4).Nodes = ( nelemx + 1)*[0:nelemy] + 1;

                                        % loop over all elements and 
                                        % set the node connectivity 
element = 0;

for j = 1 : nelemy
    for i = 1 : nelemx

        sw = i     + ( j - 1 ) * ( nelemx + 1 );
        se = i + 1 + ( j - 1 ) * ( nelemx + 1 );
        nw = i     +   j       * ( nelemx + 1 );
        ne = i + 1 +   j       * ( nelemx + 1 );

       
        if ( nnpe == 4)
            element = element + 1;

            connectivity(1,element) = sw;
            connectivity(2,element) = se;
            connectivity(3,element) = ne;
            connectivity(4,element) = nw;
            
        elseif (nnpe ==3)        % Recall that each rectangular is 
                                        % split in two triangles
            element = element + 1;

            connectivity(1,element) = sw;
            connectivity(2,element) = se;
            connectivity(3,element) = nw;

            element = element + 1;

            connectivity(1,element) = ne;
            connectivity(2,element) = nw;
            connectivity(3,element) = se;
        else
            error('This element type is not implemented\n');
        end

    end
end
connectivity=connectivity';
                                        % Set up the boundary and surface 
                                        % indicators for elems - find 
                                        % elements that have sides on the
                                        % boundaries

switch nnpe
    case 4
                                        % The first boundary (bottom)
    BoundaryElems(1).Elems = [1:nelemx];
    BoundaryElems(1).SurfaceIndicator = -2*ones(nelemx,1);

                                        % The second boundary (right)
    BoundaryElems(2).Elems = [1:nelemy]*nelemx;
    BoundaryElems(2).SurfaceIndicator = ones(nelemy,1);

                                        % The third boundary (top)
    BoundaryElems(3).Elems =  nelemx * (nelemy-1)  + [1 : nelemx];
    BoundaryElems(3).SurfaceIndicator = 2*ones(nelemx,1);

                                        % The fourth boundary (left)
    BoundaryElems(4).Elems =  nelemx *[0:(nelemy-1)] + 1;
    BoundaryElems(4).SurfaceIndicator = -1*ones(nelemy,1);

    case 3
                                        % The first boundary (bottom)
    BoundaryElems(1).Elems = [1:2:2*nelemx];
    BoundaryElems(1).SurfaceIndicator = -2*ones(nelemx,1);

                                        % The second boundary (right)
    BoundaryElems(2).Elems = [1:nelemy]*2*nelemx;
    BoundaryElems(2).SurfaceIndicator = -1*ones(nelemy,1);

                                        % The third boundary (top)
    BoundaryElems(3).Elems =  2*nelemx * (nelemy-1)+1  + [1 : 2:2*nelemx];
    BoundaryElems(3).SurfaceIndicator = -2*ones(nelemx,1);

                                        % The fourth boundary (left)
    BoundaryElems(4).Elems =  2*nelemx *[0:(nelemy-1)] + 1;
    BoundaryElems(4).SurfaceIndicator = -1*ones(nelemy,1);

    otherwise
    error('This element type is not implemented\n');
end

BoundaryStruct.elements=BoundaryElems;
BoundaryStruct.nodes=BoundaryNodes;




