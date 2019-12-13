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
    
    %      case 4  % 3x3 integration rule for Adini-Clough element
    %          % corresponding integration points and weights
    %         % 3-point quadrature rule
    %         gp = [-sqrt(3/5),-sqrt(3/5); 0,-sqrt(3/5); sqrt(3/5),-sqrt(3/5);-sqrt(3/5),0;
    %                0, 0; sqrt(3/5),0;-sqrt(3/5),sqrt(3/5); 0,sqrt(3/5);sqrt(3/5),sqrt(3/5)];
    %         w = [25/81 ;40/81 ;25/81 ; 40/81 ; 64/81 ; 40/81 ;  25/81;  40/81; 25/81];
    case 4   % 5x5 integration rule for Adini-Clough element
        % corresponding integration points and weights
        % 5-point quadrature rule
        qp(1) = - 0.906179845938663992797626878299;
        qp(2) = - 0.538469310105683091036314420700;
        qp(3) =   0.0;
        qp(4) =   0.538469310105683091036314420700;
        qp(5) =   0.906179845938663992797626878299;
        gp = [qp(1),qp(1); qp(2),qp(1);qp(3),qp(1);qp(4),qp(1);qp(5),qp(1);
            qp(1),qp(2); qp(2),qp(2);qp(3),qp(2);qp(4),qp(2);qp(5),qp(2);
            qp(1),qp(3); qp(2),qp(3);qp(3),qp(3);qp(4),qp(3);qp(5),qp(3);
            qp(1),qp(4); qp(2),qp(4);qp(3),qp(4);qp(4),qp(4);qp(5),qp(4);
            qp(1),qp(5); qp(2),qp(5);qp(3),qp(5);qp(4),qp(5);qp(5),qp(5)];
        we(1) = 0.236926885056189087514264040720;
        we(2) = 0.478628670499366468041291514836;
        we(3) = 0.568888888888888888888888888889;
        we(4) = 0.478628670499366468041291514836;
        we(5) = 0.236926885056189087514264040720;
        w = [we(1)*we(1); we(2)*we(1);we(3)*we(1);we(4)*we(1);we(5)*we(1);
            we(1)*we(2); we(2)*we(2);we(3)*we(2);we(4)*we(2);we(5)*we(2);
            we(1)*we(3); we(2)*we(3);we(3)*we(3);we(4)*we(3);we(5)*we(3);
            we(1)*we(4); we(2)*we(4);we(3)*we(4);we(4)*we(4);we(5)*we(4);
            we(1)*we(5); we(2)*we(5);we(3)*we(5);we(4)*we(5);we(5)*we(5)];     
    otherwise
        error('This element type is not implemented\n');
end


