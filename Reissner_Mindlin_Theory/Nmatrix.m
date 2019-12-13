% *************************************************************************
% Name: Nmatrix.m
% Notes: N=Nmatrix(QP)
%     This returns the matrix of shape functions N calculated at the local
%     coordinates XI,ETA for 2D elements.
%     last edit: 29 April 2015 H. Ritz
% Project Updates: Pulled from LinElast code base from MAE-5700 course,
%     updates were made to implement Kirchhoff Theory Plate Bending
% Update Authors: Dalton and Sairam
% *************************************************************************

function N=Nmatrix(LocPos)
% get xi and eta
xi=LocPos(1); eta=LocPos(2);

% shape functions depend on the number of nodes per element
% Q4 elements
N1 = (1-xi)*(1-eta)/4;
N2 = (1+xi)*(1-eta)/4;
N3 = (1+xi)*(1+eta)/4;
N4 = (1-xi)*(1+eta)/4;

% Calculate the 4 basis functions at (xi,eta). N is a row vector
N = [N1, N2, N3, N4];  
end


