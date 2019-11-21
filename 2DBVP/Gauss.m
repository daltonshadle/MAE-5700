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

function [gp,w] = Gauss(nnpe)

%  Provides Gauss points and weights for different elements
%  Additional info needs to be added for other than the Q4 and T3 elements.

switch nnpe
    case 4 % 2x2 integration rule for Q4 element
        pt=1/sqrt(3);
        gp = [-pt,-pt; -pt,pt; pt,-pt; pt,pt];
        w  = [1,1,1,1];
        
    case 3 % The T3 element
        pt=1/6;
        gp =[pt,pt; 4*pt,pt; pt,4*pt];  
        w = [1/6 1/6 1/6];      
    otherwise
        error('This element type is not implemented\n');
end


