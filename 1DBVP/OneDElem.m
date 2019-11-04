function [ke,fe] = OneDElem(elmID,qp,w,meshStruct,boundStruct)
% [localstiffnessmatrix, localforcevector] = OneDElem(elementnumber, QuadPoints,Weights,meshStruct)
% generate the local stiffness matrix and local force vector from
% body forces for use with ONEDBVP code.
% last edit: 20 October 2015 H. Ritz

% unpack necessary input
elCon  =meshStruct.elCon;
nCoords=meshStruct.nCoords;
numDOF =meshStruct.numDOF;
nnpe   =meshStruct.nnpe;

appForces=boundStruct.appForces;

ke=zeros(numDOF*nnpe); % initialize elemental stiffness matrix
fe=zeros(numDOF*nnpe,1); % initialize elemental force vector

% get nodal coordinates for this element
xvec=nCoords(elCon(elmID,:))'; % this is a column vector of the
                               % physical coordinates of the
                               % nodes of this element

for iqp=1:length(qp) % loop over quadrature points
    % for each quadrature point ...
    NofXi=Nmatrix(qp(iqp),nnpe); % row vector of shape functions
    % at this QP in parent domain
    dNdXi=Bmatrix(qp(iqp),nnpe); % row vector of shape function derivatives
    % at this QP in parent domain
    JofXi=dNdXi*xvec; % jacobian at this QP (scalar)
    B=dNdXi/JofXi; % row vector of dNdX at ths QP, needed for integration
    
    X=NofXi*xvec; % physical coordinate of this QP (scalar)
    
    % add the weighted contributions of this QP to the elemental stiffness
    % matrix and elemental body force vector
    ke=ke+(PP(X)*B'*B)*w(iqp)*JofXi;
    fe=fe+(NofXi'*FF(X))*w(iqp)*JofXi;
end

% check if there are any point forces inside this element
for ipf=1:size(appForces,1) % for each point force
    xp=appForces(ipf,1);    % location of point force
    val=appForces(ipf,2);   % value of point force
    if ((xp>min(xvec))&&(xp<=max(xvec))) % if the x position of the force is within this element
        xip=2*(xp-xvec(1))/(xvec(end)-xvec(1))-1; % find the local coordinate of xp
        Np=Nmatrix(xip,nnpe);    % evaluate the shape functions here
        fe=fe+Np'*val;      % spread the point force to the nodes and add to the force vector
    end
end
