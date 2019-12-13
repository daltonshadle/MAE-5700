function ddN=ddNmatrix(LocPos)
% N=Nmatrix(QP)
% This returns the matrix of shape functions N calculated at the local
% coordinates XI,ETA for 2D elements.
% last edit: 29 April 2015 H. Ritz
% Calculate the derivatives of the basis functions with respect to xi and eta. dN is a 2-row vector
xi=LocPos(1); eta=LocPos(2);
% derivatives of shape function w.r.t to xi
ddN_by_xi_1 = -3*xi*(eta-1)/4;
ddN_by_xi_2 = -(6*xi-2)*(eta-1)/8;
ddN_by_xi_3 = 0;
ddN_by_xi_4 = 3*xi*(eta-1)/4;
ddN_by_xi_5 = -(6*xi+2)*(eta-1)/8;
ddN_by_xi_6 = 0;
ddN_by_xi_7 = -3*xi*(eta+1)/4;
ddN_by_xi_8 = (6*xi+2)*(eta+1)/8;
ddN_by_xi_9 = 0;
ddN_by_xi_10 = 3*xi*(eta+1)/4;
ddN_by_xi_11 = (6*xi-2)*(eta+1)/8;
ddN_by_xi_12 = 0;
% derivatives of shape function w.r.t to eta
ddN_by_eta_1 = -3*eta*(xi-1)/4;
ddN_by_eta_2 = 0;
ddN_by_eta_3 = -(6*eta-2)*(xi-1)/8;
ddN_by_eta_4 = 3*eta*(xi+1)/4;
ddN_by_eta_5 = 0;
ddN_by_eta_6 = (6*eta-2)*(xi+1)/8;
ddN_by_eta_7 = -3*eta*(xi+1)/4;
ddN_by_eta_8 = 0;
ddN_by_eta_9 = (6*eta+2)*(xi+1)/8;
ddN_by_eta_10 = 3*eta*(xi-1)/4;
ddN_by_eta_11 = 0;
ddN_by_eta_12 = -(6*eta+2)*(xi-1)/8;
% derivatives of shape function w.r.t to eta and xi
ddN_by_xi_eta_1 = -(3/8)*eta^2 - (3/8)*xi^2 + 1/2;
ddN_by_xi_eta_2 = -(3/8)*xi^2 + (1/4)*xi + 1/8;
ddN_by_xi_eta_3 =  -3*eta^2/8 + eta/4 + 1/8;
ddN_by_xi_eta_4 =  (3/8)*eta^2 + (3/8)*xi^2 - 1/2;
ddN_by_xi_eta_5 =   -(3/8)*xi^2 - (1/4)*xi + 1/8;
ddN_by_xi_eta_6 =   3*eta^2/8 - eta/4 - 1/8;
ddN_by_xi_eta_7 =   -(3/8)*eta^2 - (3/8)*xi^2 + 1/2;
ddN_by_xi_eta_8 =  (3/8)*xi^2 + (1/4)*xi - 1/8;
ddN_by_xi_eta_9 =  3*eta^2/8 + eta/4 - 1/8;
ddN_by_xi_eta_10 =  (3/8)*eta^2 + (3/8)*xi^2 - 1/2;
ddN_by_xi_eta_11 = (3/8)*xi^2 -xi/4 -1/8;
ddN_by_xi_eta_12 = (-3/8)*eta^2 -(1/4)*eta + (1/8);
% derivatives placed together
ddN = [ddN_by_xi_1, ddN_by_xi_2, ddN_by_xi_3, ddN_by_xi_4, ddN_by_xi_5, ddN_by_xi_6, ...
    ddN_by_xi_7, ddN_by_xi_8, ddN_by_xi_9, ddN_by_xi_10, ddN_by_xi_11, ddN_by_xi_12;
    ddN_by_eta_1, ddN_by_eta_2, ddN_by_eta_3, ddN_by_eta_4, ddN_by_eta_5, ddN_by_eta_6, ...
    ddN_by_eta_7, ddN_by_eta_8, ddN_by_eta_9, ddN_by_eta_10, ddN_by_eta_11, ddN_by_eta_12;    
    ddN_by_xi_eta_1, ddN_by_xi_eta_2, ddN_by_xi_eta_3, ddN_by_xi_eta_4, ddN_by_xi_eta_5, ddN_by_xi_eta_6, ...
    ddN_by_xi_eta_7, ddN_by_xi_eta_8, ddN_by_xi_eta_9, ddN_by_xi_eta_10, ddN_by_xi_eta_11, ddN_by_xi_eta_12];
end
