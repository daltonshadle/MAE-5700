% ------------------------------------------------------------------------|
%                                                                         |
% MAE4700-5700, Finite Element Analysis for Mechanical & Aerospace Design |
%                                                                         |
% Copyright: Cornell University (this software should not be used without |
% written permission)                                                     |
%                                                                         |
% Authors: N. Zabaras (zabaras@cornell.edu) & Xiang Ma (xm25@cornell.edu) |
% Last update: 17 Nov 2015 H. Ritz ; Y. Xu                                |
% ------------------------------------------------------------------------|
%
function  globalSystem=ApplyNaturalBC(i,boundStruct,meshStruct,globalSystem)
% Computed the element load from natural BCs.

% unpack necessary variables
sideInd = boundStruct.SurfNat(i,1);
qt = boundStruct.SurfNat(i,2); % tangetial traction
qn = boundStruct.SurfNat(i,3); % normal traction
BElems = boundStruct.elements(sideInd).Elems;
SurfID = boundStruct.elements(sideInd).SurfaceIndicator;

nnpe=meshStruct.nnpe;
numDOF=meshStruct.numDOF;
nCoords=meshStruct.nCoords;
elCon=meshStruct.elCon;
gatherMat=meshStruct.gatherMat;
% K = globalSystem.K;
F = globalSystem.F;
ldof = nnpe*numDOF;       % Degrees of freedom per element

Num = length(BElems); %  Extract the number of boundary elements

pt=1/sqrt(3);          % Set up Gauss quadrature for boundary integral
gp = [-pt,pt];         % 2 points for one dimension - this will need to
w  = [1,1];            % be generalized for other situations.

for elmID = 1 : Num
    
%     Ke = zeros(ldof);     % initialize element stiffness matrix
    fe = zeros(ldof,1);   % initialize element load vector
    
    glb = elCon(BElems(elmID),:);   % Global number of the element nodes
    
    coord = nCoords(glb,:);          % Global nodal coordinates for the element
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
        
        if ( nnpe == 4)
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
            
        elseif  ( nnpe == 3)
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
        
        felm.nx = Jacc(2)/detJ;    % Calculate the direct cosines
        felm.ny = -Jacc(1)/detJ;

        % Extract the N matrix
        if (nnpe == 4)
            N = [ felm.N(1),   0,    felm.N(2),   0,     felm.N(3),     0,  felm.N(4),  0;
                0 ,  felm.N(1),    0   ,  felm.N(2),   0,    felm.N(3),  0 ,  felm.N(4)];
        else
            N = [ felm.N(1),   0,    felm.N(2),   0,     felm.N(3),     0;
                0 ,  felm.N(1),    0   ,  felm.N(2),   0,    felm.N(3)];
            
        end
        
        % Transform to x and y direction using the following formula
        % where nx and ny are the direction cosines at the current boundary
        q = [felm.nx -felm.ny;felm.ny  felm.nx] * [qn;qt];
        
        % Assemble into force vector
        fe = fe +  N' * q * felm.detJxW;
        
        
    end
    clear felm;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Assemble to global
    glbROW = gatherMat(BElems(elmID),:); % global row index
    F(glbROW)=F(glbROW)+fe;

end
% globalSystem.K = K;
globalSystem.F = F;