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
dN_by_xi_1 = (eta-1)/4;
dN_by_xi_2 = (1-eta)/4;
dN_by_xi_3 = (1+eta)/4;
dN_by_xi_4 = -(1+eta)/4;

% derivatives of shape function w.r.t to eta
dN_by_eta_1 = (xi-1)/4; 
dN_by_eta_2 = -(xi+1)/4;
dN_by_eta_3 = (1+xi)/4; 
dN_by_eta_4 = (1-xi)/4;

% derivatives placed together, dN (2x4 matrix)
dN = [dN_by_xi_1, dN_by_xi_2, dN_by_xi_3, dN_by_xi_4;
      dN_by_eta_1, dN_by_eta_2, dN_by_eta_3, dN_by_eta_4];
end
