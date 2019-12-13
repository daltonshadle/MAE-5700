% *************************************************************************
% Name: dNmatrix.m
% Notes: dN=Nmatrix(QP)
%     This returns the matrix of the derivatives of shape functions dN 
%     calculated at the local coordinates XI,ETA for 2D elements. This is a
%     new file for this project.
% Authors: Dalton and Sairam
% *************************************************************************

function dN=dNmatrix(LocPos)
% Calculate the derivatives of the shape functions with respect to xi and eta.

% get xi and eta
xi=LocPos(1); eta=LocPos(2);

% derivatives of shape function w.r.t to xi
dN_by_xi_1 = (-1/8)*(eta-1)*(eta^2 + eta + 3*xi^2 -3);
dN_by_xi_2 = (1/8)*(eta-1)*(-3*xi^2 + 2*xi + 1);
dN_by_xi_3 = (-1/8)*(eta-1)^2*(eta+1);
dN_by_xi_4 = (1/8)*(eta-1)*(eta^2 + eta + 3*xi^2 -3);
dN_by_xi_5 = (-1/8)*(eta-1)*(3*xi^2 + 2*xi -1); 
dN_by_xi_6 = (1/8)*(eta-1)^2*(eta+1);
dN_by_xi_7 = (1/8)*(eta+1)*(-eta^2 + eta -3*xi^2 +3); 
dN_by_xi_8 = (1/8)*(eta+1)*(3*xi^2 + 2*xi -1);
dN_by_xi_9 = (1/8)*(eta-1)*(eta+1)^2; 
dN_by_xi_10 = (-1/8)*(eta+1)*(-eta^2 + eta - 3*xi^2 + 3);
dN_by_xi_11 = (-1/8)*(eta+1)*(-3*xi^2 + 2*xi +1);
dN_by_xi_12 = (-1/8)*(eta-1)*(eta+1)^2;

% derivatives of shape function w.r.t to eta
dN_by_eta_1 = (-1/8)*(xi-1)*(3*eta^2 + xi^2 + xi -3); 
dN_by_eta_2 = (-1/8)*(xi-1)^2*(xi+1);
dN_by_eta_3 = (1/8)*(xi-1)*(-3*eta^2 + 2*eta +1); 
dN_by_eta_4 = (-1/8)*(xi+1)*(-3*eta^2 -xi^2 +xi+3);
dN_by_eta_5 = (-1/8)*(xi-1)*(xi+1)^2; 
dN_by_eta_6 = (-1/8)*(xi+1)*(-3*eta^2 + 2*eta+1);
dN_by_eta_7 = (1/8)*(xi+1)*(-3*eta^2 -xi^2 +xi +3); 
dN_by_eta_8 = (1/8)*(xi-1)*(xi+1)^2;
dN_by_eta_9 = (1/8)*(xi+1)*(3*eta^2 + 2*eta -1); 
dN_by_eta_10 = (1/8)*(xi-1)*(3*eta^2 + xi^2 + xi -3);
dN_by_eta_11 = (1/8)*(xi-1)^2*(xi+1); 
dN_by_eta_12 = (-1/8)*(xi-1)*(3*eta^2 + 2*eta -1);

% derivatives placed together, dN (2x12 matrix)
dN = [dN_by_xi_1, dN_by_xi_2, dN_by_xi_3, dN_by_xi_4, dN_by_xi_5, dN_by_xi_6, ...
    dN_by_xi_7, dN_by_xi_8, dN_by_xi_9, dN_by_xi_10, dN_by_xi_11, dN_by_xi_12;
    dN_by_eta_1, dN_by_eta_2, dN_by_eta_3, dN_by_eta_4, dN_by_eta_5, dN_by_eta_6, ...
    dN_by_eta_7, dN_by_eta_8, dN_by_eta_9, dN_by_eta_10, dN_by_eta_11, dN_by_eta_12];
end
