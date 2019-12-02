function N=Nmatrix(LocPos)
% N=Nmatrix(QP)
% This returns the matrix of shape functions N calculated at the local
% coordinates XI,ETA for 2D elements.
% last edit: 29 April 2015 H. Ritz

xi=LocPos(1); eta=LocPos(2);
% shape functions depend on the number of nodes per element
% Q4 elements elements
N1 = (-1/8)*(eta-1)*(xi-1)*(eta^2 + eta + xi^2 + xi -2);
N2 = (-1/8)*(eta-1)*(xi-1)^2*(xi+1);
N3 = (-1/8)*(eta-1)^2*(eta+1)*(xi-1);
N4 = (1/8)*(eta-1)*(xi+1)*(eta^2 + eta + xi^2 -xi-2);
N5 = (-1/8)*(eta-1)*(xi-1)*(xi+1)^2;
N6 = (1/8)*(eta-1)^2*(eta+1)*(xi+1);
N7 = (1/8)*(eta+1)*(xi+1)*(-eta^2 +  eta - xi^2 + xi +2);
N8 = (1/8)*(eta+1)*(xi-1)*(xi+1)^2;
N9 = (1/8)*(eta-1)*(eta+1)^2*(xi+1);
N10 = (1/8)*(eta+1)*(xi-1)*(eta^2 - eta + xi^2 + xi -2);
N11 = (1/8)*(eta+1)*(xi-1)^2*(xi+1);
N12 = (-1/8)*(eta-1)*(eta+1)^2*(xi-1);
N = [N1, N2, N3, N4, N5, N6, N7, N8, N9, N10, N11, N12];  % Calculate the 12 basis functions
% at (xi,eta). N is a row
% vector
end


