function [ke,fe] = TwoDElem(MeshStruct,elmID,qp,w)
% [localstiffnessmatrix, localforcevector] = TwoDElem(elementnumber,QuadPoints,Weights)
% generate the local stiffness matrix and local force vector from
% body forces for use with TWODBVP code.
% last edit: 30 April 2015 H. Ritz

% unpack necessary information
nnpe=MeshStruct.nnpe;
numDOF=MeshStruct.numDOF;
nCoords=MeshStruct.nCoords;
elCon=MeshStruct.elCon;

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
    B=inv(JofXiEta)*dNdXiEta;         % 2xnnpe array of dNdX at ths QP
    
    XY=NofXiEta*[xvec,yvec]; % physical coordinate of this QP (1x2)
    
    % add the weighted contributions of this QP to the elemental stiffness
    % matrix and elemental body force vector
    ke=ke+(B'*DD(XY)*B)*w(iqp)*det(JofXiEta);
    fe=fe+(FF(XY)*NofXiEta')*w(iqp)*det(JofXiEta);
end