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

function value = FF ( x )

% It returns the rhs f(x,y) of the PDE at point x=(x(1),x(2)): 
%  -grad . ([k]gradu) =f(x,y) or for [k] diagonal matrix, 
% - d/dx(kx*du/dx) - d/dy(ky*du/dy)  =  f(x,y)
%
% This is problem specific. 

% value = 0;

% For Torsion
G = 77.2e+9; % N/m^2
alpha = 0.3;
value = 2 * G * alpha;

% For Heat Transfer
value=0;