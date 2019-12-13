% *************************************************************************
% Name: TwoDElem.m
% Notes: [localstiffnessmatrix, localforcevector] = 
%     TwoDElem(meshStruct,elementnumber,QuadPoints,Weights)
%     generate the local stiffness matrix and local force vector from
%     body forces for use with LINELAST code.
%     last edit: 1 May 2015 H. Ritz
% Project Updates: Pulled from LinElast code base from MAE-5700 course,
%     updates were made to implement Kirchhoff Theory Plate Bending
% Update Authors: Dalton and Sairam
% *************************************************************************

function [ke,fe] = TwoDElem(meshStruct,elmID,qp,w)
%% unpack necessary information
nnpe=meshStruct.nnpe;
numDOF=meshStruct.numDOF;
nCoords=meshStruct.nCoords;
elCon=meshStruct.elCon;
Db=meshStruct.Material.Db;
Ds=meshStruct.Material.Ds;
sc=meshStruct.Material.sc;
thickness = meshStruct.thickness;

ke=zeros(numDOF*nnpe);   % initialize elemental stiffness matrix
fe=zeros(numDOF*nnpe,1); % initialize elemental force vector

%% get nodal coordinates for this element
xvec=nCoords(elCon(elmID,:),1); % these are column vectors of the
yvec=nCoords(elCon(elmID,:),2); % physical coordinates of the
                                % nodes of this element
                                

% iterate over each quadrature point
for iqp=1:size(qp,1)
    % for each quadrature point ...
    NofXiEta=Nmatrix(qp(iqp,:)); % row vector of shape functions
                                      % at this QP in parent domain (1x4)
    dNdXiEta=dNmatrix(qp(iqp,:));% 2xnnpe array of shape func derivs
                                      % at this QP in parent domain (2x4)
    JofXiEta=dNdXiEta*[xvec,yvec];    % jacobian at this QP (2x2 matrix)
    dNdXY=inv(JofXiEta)*dNdXiEta;     % 2xnnpe array of dNdX at ths QP
    
    XY=NofXiEta*[xvec,yvec]; % physical coordinate of this QP (1x2)
    
    % now create the N and Bb (bending) and Bs (shear) matrices
    N=[]; Bb=[]; Bs=[];
    for np=1:nnpe
        N=[N,[NofXiEta(np), 0, 0; 0, NofXiEta(np), 0; 0, 0, NofXiEta(np)]];
        Bb=[Bb,[0, dNdXY(1,np), 0; 0, 0, dNdXY(2,np); 0, dNdXY(2,np), dNdXY(1,np)]];
        Bs=[Bs,[dNdXY(1,np), -NofXiEta(np), 0; dNdXY(2,np), 0, -NofXiEta(np)]];
    end
    
    % add the weighted contributions of this QP to the elemental stiffness
    % matrix and elemental body force vector
    ke=ke+thickness^3/12*(Bb'*Db*Bb)*w(iqp)*det(JofXiEta)+sc*thickness*(Bs'*Ds*Bs)*w(iqp)*det(JofXiEta);
    fe=fe+(N'*BodyForce(XY))*w(iqp)*det(JofXiEta);
end