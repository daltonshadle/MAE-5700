% ------------------------------------------------------------------------|
%                                                                         |
% MAE4700-5700, Finite Element Analysis for Mechanical & Aerospace Design |
%                                                                         |
% Copyright: Cornell University (this software should not be used without |
% written permission)                                                     |
%                                                                         |
% Authors: N. Zabaras (zabaras@cornell.edu) & Xiang Ma (xm25@cornell.edu) |
%                                                                         |
% ------------------------------------------------------------------------|
%

function [gp,w] = OneDGauss(nnpe)

%  Provides Gauss points and weights for different elements
%  Additional info needs to be added for other than the Q4 and T3 elements.

switch nnpe
    case 4   % 5x5 integration rule for Adini-Clough element
        % corresponding integration points and weights
        % 5-point quadrature rule
        qp(1) = - 0.906179845938663992797626878299;
        qp(2) = - 0.538469310105683091036314420700;
        qp(3) =   0.0;
        qp(4) =   0.538469310105683091036314420700;
        qp(5) =   0.906179845938663992797626878299;
        we(1) = 0.236926885056189087514264040720;
        we(2) = 0.478628670499366468041291514836;
        we(3) = 0.568888888888888888888888888889;
        we(4) = 0.478628670499366468041291514836;
        we(5) = 0.236926885056189087514264040720;
        
        gp = qp;
        w = we;
    otherwise
        error('This element type is not implemented\n');
end


