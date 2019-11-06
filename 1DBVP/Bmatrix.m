function B=Bmatrix(xi,nnpe)
% B=Bmatrix(XI,nnpe)
% This returns the matrix of shape function derivatives calculated at the 
% local coordinate XI for 1D elements.
% last edit: 6 August 2015 H. Ritz

switch nnpe % shape functions depend on the number of nodes per element
    case 2 % linear elements
        B=(1/2)*[-1, 1];
    case 3 % quadratic elements
        B=[xi-(1/2), -2*xi, xi+(1/2)];
    case 4 % cubic elements
        B = [(9/16)*(-3*xi^2 + 2*xi + 1/9), (27/16)*(3*xi^2 - 2*xi/3 - 1), ...
            (27/16)*(-3*xi^2 - 2*xi/3 + 1), (9/16)*(3*xi^2 + 2*xi - 1/9)];
end
