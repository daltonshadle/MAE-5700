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
                                
x_elemlength =  abs(xvec(2) - xvec(1))/2;                           
y_elemlength =  abs(yvec(3) - yvec(2))/2;         
x_elem_center = mean(xvec);
y_elem_center = mean(yvec);

for iqp=1:size(qp,1) % loop over quadrature points
    % for each quadrature point ...
    NofXiEta=Nmatrix(qp(iqp,:));  % 1x12 matrix,  row vector of shape functions
                                      % at this QP in parent domain
    dNdXiEta=dNmatrix(qp(iqp,:)); % 2x12 matrix, of first derivs of shape func 
                                      % at this QP in parent domain                                  
    ddNdXiEta=ddNmatrix(qp(iqp,:)); % 3x12 matrix, of second derivs of shape func 
                                      % at this QP in parent domain                                    
                                          
    JofXiEta= [x_elemlength, 0; 0, y_elemlength]; % jacobian at this QP (2x2 matrix)

    XY= [x_elemlength*qp(iqp,1)+x_elem_center, y_elemlength*qp(iqp,2)+y_elem_center];   % physical coordinate of this QP (1x2)
    
    % now create the N and B matrices
    
    N = [NofXiEta; (1/x_elemlength)*dNdXiEta(1,:); (1/y_elemlength)*dNdXiEta(2,:)];
    
    B =  [(1/x_elemlength)^2*ddNdXiEta(1,:); (1/y_elemlength)^2*ddNdXiEta(2,:); (2/x_elemlength)*(1/y_elemlength)*ddNdXiEta(3,:)];
    
    % add the weighted contributions of this QP to the elemental stiffness
    % matrix and elemental body force vector
    ke=ke+(B'*D*B)*w(iqp)*det(JofXiEta);
    fe=fe+(N'*BodyForce(XY))*w(iqp)*det(JofXiEta);
end