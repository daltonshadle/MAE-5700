% ------------------------------------------------------------------------|
%                                                                         |
% MAE4700-5700, Finite Element Analysis for Mechanical & Aerospace Design |
%                                                                         |
% Copyright: Cornell University (this software should not be used without |
% written permission)                                                     |
%                                                                         |
% Authors: N. Zabaras (zabaras@cornell.edu) & Xiang Ma (xm25@cornell.edu) |
% Last update: 1 May 2015 H. Ritz                                         |                                  |
% ------------------------------------------------------------------------|
%

function value = BodyForce ( x )

% It returns the body force vector at every location.
%
% This is problem specific. 

 mua=1e7; % mass per unit area
% value = [0;-9.81*mua]; % use this for gravity in the -y direction

value=[-9.81*mua; 0; 0]; % use this for no body force
