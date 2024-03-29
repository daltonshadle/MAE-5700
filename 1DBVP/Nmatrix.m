function N=Nmatrix(xi,nnpe)
% N=Nmatrix(XI,nnpe)
% This returns the matrix of shape functions N calculated at the local
% coordinate XI for 1D elements.
% last edit: 6 August 2015 H. Ritz

switch nnpe % shape functions depend on the number of nodes per element
    case 2 % linear elements
        N=(1/2)*[(1-xi), (xi+1)];
    case 3 % quadratic elements
        N=[(1/2)*xi*(xi-1), -(1+xi)*(-1+xi), (1/2)*xi*(xi+1)];
    case 4 % cubic elements
        N = [(9/16)*(1-xi)*(xi^2- 1/9), (27/16)*(xi^2-1)*(xi-1/3), ...
            (27/16)*(1-xi^2)*(xi+1/3), (9/16)*(1+xi)*(xi^2- 1/9)];
end
