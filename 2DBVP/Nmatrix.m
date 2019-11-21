function N=Nmatrix(LocPos,nnpe)
% N=Nmatrix(QP)
% This returns the matrix of shape functions N calculated at the local
% coordinates XI,ETA for 2D elements.
% last edit: 29 April 2015 H. Ritz

xi=LocPos(1); eta=LocPos(2);
switch nnpe % shape functions depend on the number of nodes per element
    case 4 % Q4 elements elements
    N = 1/4*[(1 - xi).*(1-eta),(1+xi).*(1-eta),...
              (1+xi).*(1+eta), (1-xi).*(1+eta)];  % Calculate the 4 basis functions 
                                                  % at (xi,eta). N is a row
                                                  % vector
    case 3 % T3 elements elements
    N = [xi, eta, 1-xi-eta];                      % Calculate the 3 basis functions 
                                                  % at (xi,eta). N is a row
                                                  % vector
end
