function [ke,fe] = TwoDElem(meshStruct,elmID,qp,w)
% [localstiffnessmatrix, localforcevector] = TwoDElem(meshStruct,elementnumber,QuadPoints,Weights)
% generate the local stiffness matrix and local force vector from
% body forces for use with LINELAST code.
% last edit: 1 May 2015 H. Ritz

% unpack necessary information
nnpe=meshStruct.nnpe;
numDOF=meshStruct.numDOF;
nCoords=meshStruct.nCoords;
elCon=meshStruct.elCon;
D=meshStruct.Material.D;

ke=zeros(numDOF*nnpe);   % initialize elemental stiffness matrix
fe=zeros(numDOF*nnpe,1); % initialize elemental force vector

% get nodal coordinates for this element
xvec=nCoords(elCon(elmID,:),1); % these are column vectors of the
yvec=nCoords(elCon(elmID,:),2); % physical coordinates of the
                                % nodes of this element
                               

for iqp=1:size(qp,1) % loop over quadrature points
    % for each quadrature point ...
    NofXiEta=Nmatrix(qp(iqp,:),nnpe); % row vector of shape functions
                                      % at this QP in parent domain
    dNdXiEta=dNmatrix(qp(iqp,:),nnpe);% 2xnnpe array of shape func derivs
                                      % at this QP in parent domain
    JofXiEta=dNdXiEta*[xvec,yvec];    % jacobian at this QP (2x2 matrix)
    dNdXY=inv(JofXiEta)*dNdXiEta;     % 2xnnpe array of dNdX at ths QP
    
    XY=NofXiEta*[xvec,yvec]; % physical coordinate of this QP (1x2)
    
    % now create the N and B matrices
    N=[]; B=[];
    for np=1:nnpe
        N=[N,[NofXiEta(np), 0; 0, NofXiEta(np)]];
        B=[B,[dNdXY(1,np), 0; 0, dNdXY(2,np); dNdXY(2,np), dNdXY(1,np)]];
    end
    % add the weighted contributions of this QP to the elemental stiffness
    % matrix and elemental body force vector
    ke=ke+(B'*D*B)*w(iqp)*det(JofXiEta);
    fe=fe+(N'*BodyForce(XY))*w(iqp)*det(JofXiEta);
end