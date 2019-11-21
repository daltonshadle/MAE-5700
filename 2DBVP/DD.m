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
function D = DD(XY)

% Sets up the `conductivity' matrix (2x2 for 2D) 
% at point XY=(x(1),x(2)).
%
% This is problem specific (constitutive) information relating 
% the flux and the gradient of the main field variable.

% D = 5*[1  0;
%      0  1];
D = [1 0;
     0 1;];