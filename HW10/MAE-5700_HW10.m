%**************************************************************************
% Name: MAE-5700_HW10.m
% Author: Dalton and Sairam
% Date: Nov. 15, 2019
% Notes: Summarized changes to HW10 Question 2 in HW7Problem2.m
%**************************************************************************

% From HW7Problem2.m ******************************************************
% ...
for row=1:size(XI,1)
    for row=1:size(XI,1)
        % ...
        x(row,col)= N * XY(:,1); % use the shape functions to find the global 
        y(row,col)= N * XY(:,2); % coordinates of this \xi, \eta
    end
end

% ...
function JacCalcFun(XY)
switch nnpe
    % ...
    case 8
        % Shape function for Q8 condition, results in [8x1] matrix
        N=[(-1/4)*(1-xi)*(1-eta)*(1+xi+eta);
           (-1/4)*(1+xi)*(1-eta)*(1-xi+eta);
           (-1/4)*(1+xi)*(1+eta)*(1-xi-eta);
           (-1/4)*(1-xi)*(1+eta)*(1+xi-eta);
           (1/2)*(1-xi^2)*(1-eta);
           (1/2)*(1+xi)*(1-eta^2);
           (1/2)*(1-xi^2)*(1+eta);
           (1/2)*(1-xi)*(1-eta^2)]';
end

% ...
function GN=GNmat(xi, eta, nnpe)
switch nnpe
    case 4
        % Shape function gradient for Q4 condition, results in [2x4] matrix
        GN=(1/4)*[(eta-1), (1-eta), (eta+1), -(eta+1);
                  (xi-1),  (1-xi),  (xi+1),  -(xi+1)];
    % ...
end

% ...
function Jac=JacMat(xi, eta, XY)
    GN=GNmat(xi,eta,size(XY,1));
    Jac=GN * XY;
