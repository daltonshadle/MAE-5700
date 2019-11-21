% ------------------------------------------------------------------------|
%                                                                         |
% MAE4700-5700, Finite Element Analysis for Mechanical & Aerospace Design |
%                                                                         |
% Copyright: Cornell University (this software should not be used without |
% written permission)                                                     |
%                                                                         |
% Authors: N. Zabaras (zabaras@cornell.edu) & Xiang Ma (xm25@cornell.edu) |
% Last update: 7 November 2017 H. Ritz                                      |                                  |
% ------------------------------------------------------------------------|
%
function  globalSystem=ApplyNaturalBC(i,boundStruct,meshStruct,globalSystem)
% Compute the element stiffness and load from natural BCs.
sideInd = boundStruct.SurfNat(i,1);
alpha = boundStruct.SurfNat(i,2);
beta = boundStruct.SurfNat(i,3);

% unpack necessary variables
BoundaryElems=boundStruct.elements;
nen=meshStruct.nnpe;
ndof=meshStruct.numDOF;
Elems=meshStruct.elCon';
Nodes=meshStruct.nCoords';
gatherMat=meshStruct.gatherMat;
K = globalSystem.K;
F = globalSystem.F;
ldof = nen*ndof;       % Degrees of freedom per element

%  Extract the element number on boundary ID.

BElems = BoundaryElems(sideInd).Elems;
SurfID = BoundaryElems(sideInd).SurfaceIndicator;

%  Extract the number of boundary elements

Num = length(BElems);

pt=1/sqrt(3);          % Set up Gauss quadrature for boundary integral
gp = [-pt,pt];         % 2 points for one dimension - this will need to
w  = [1,1];            % be generalized for other situations.

for elmID = 1 : Num
    
    Ke = zeros(ldof);     % initialize element stiffness matrix
    fe = zeros(ldof,1);   % initialize element load vector
    
    glb = Elems(:,BElems(elmID));   % Global number of the element nodes
    
    coord = Nodes(:,glb)';          % Global nodal coordinates for the element
    % nodes as a column vector
    
    for qp = 1:length(w)             % loop over all the gauss points on the
        % boundary segment
        
        %--------------------------------------------------------------------
        % First set up the information of the finite element for the current
        % integration point. In other words, put the value and derivative of
        % basis function, the physical coordinates of the integration points in
        % the struct fe.
        %------------------------------------------------------------------
        %
        % Different basis functions according to different sides
        
        if ( nen == 4)
            switch (SurfID(elmID))
                case -2
                    n = [(1 - gp(qp))/2, (1 + gp(qp))/2, 0, 0];
                    dn = [-1/2, 1/2, 0, 0];
                case +1
                    n = [0, (1 - gp(qp))/2, (1 + gp(qp))/2, 0];
                    dn = [0, -1/2, 1/2, 0];
                case +2
                    n = [0, 0, (1 - gp(qp))/2, (1 + gp(qp))/2];
                    dn = [0, 0, -1/2, 1/2];
                case -1
                    n = [(1 + gp(qp))/2, 0, 0, (1 - gp(qp))/2];
                    dn = [1/2, 0, 0, -1/2];
            end
            
        elseif  ( nen == 3)
            switch (SurfID(elmID))
                case -2
                    n = [(1 - gp(qp))/2, (1 + gp(qp))/2, 0];
                    dn = [-1/2, 1/2, 0];
                case +1
                    n = [0, (1 - gp(qp))/2, (1 + gp(qp))/2];
                    dn = [0, -1/2, 1/2];
                case -1
                    n = [(1 + gp(qp))/2, 0, (1 - gp(qp))/2];
                    dn = [1/2, 0, -1/2];
            end
        else
            fprintf(1, 'This element type is not implemented\n');
        end
        
        felm.N = n;           % put basis function into struct fe
        felm.x = n*coord;     % Calculate the global coordinates of the
        % current integration points
        Jacc = dn*coord;      % Calculate the components of the Jaccobian matrix
        % Now it is a vector
        detJ = norm(Jacc);    % Calculate the norm of the vector
        
        felm.detJxW = detJ * w(qp); % Calculate the integration weight
        
        if ( alpha ~=0)
            Ke = Ke - felm.N' * felm.N * alpha * felm.detJxW;
        end
        
        if ( beta ~=0)
            fe = fe + felm.N' * beta * felm.detJxW;
        end
        
        
    end
    clear felm;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nnpe=nen;
    numDOF=ndof;
    elCon=Elems';
        for Lrw = 1 : (nnpe*numDOF)
            Grw = gatherMat(BElems(elmID),Lrw); % global row index
            F(Grw) = F(Grw) + fe(Lrw);
            for Lcl = 1 : (nnpe*numDOF)
                Gcl = gatherMat(BElems(elmID),Lcl); % global column index
                K(Grw,Gcl) = K(Grw,Gcl) + Ke(Lrw,Lcl);
                
            end
        end
end
globalSystem.K = K;
globalSystem.F = F;