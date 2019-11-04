function [qp,w] = Gauss(nqp)
% [GaussPointLocation,Weight]=GAUSS(NumberOfGaussPoints)
%  Returns the local gauss quadrature point locations and associated
%  weights for 1D problems.
switch nqp
    case 1   % 1 Gauss point
        qp = 0;
        w  = 2;
        
    case 2   % 2 Gauss points
        qp(1) = - 0.577350269189625764509148780502;
        qp(2) =   0.577350269189625764509148780502;
        
        w(1) = 1.0;
        w(2) = 1.0;
        
    case 3   % 3 Gauss points
        qp(1) = - 0.774596669241483377035853079956;
        qp(2) =   0.0;
        qp(3) =   0.774596669241483377035853079956;
        
        w(1) = 5.0 / 9.0;
        w(2) = 8.0 / 9.0;
        w(3) = 5.0 / 9.0;
        
    case 4   % 4 Gauss points        
        qp(1) = - 0.861136311594052575223946488893;
        qp(2) = - 0.339981043584856264802665759103;
        qp(3) =   0.339981043584856264802665759103;
        qp(4) =   0.861136311594052575223946488893;
        
        w(1) = 0.347854845137453857373063949222;
        w(2) = 0.652145154862546142626936050778;
        w(3) = 0.652145154862546142626936050778;
        w(4) = 0.347854845137453857373063949222;
        
    case 5   % 5 Gauss points
        qp(1) = - 0.906179845938663992797626878299;
        qp(2) = - 0.538469310105683091036314420700;
        qp(3) =   0.0;
        qp(4) =   0.538469310105683091036314420700;
        qp(5) =   0.906179845938663992797626878299;
        
        w(1) = 0.236926885056189087514264040720;
        w(2) = 0.478628670499366468041291514836;
        w(3) = 0.568888888888888888888888888889;
        w(4) = 0.478628670499366468041291514836;
        w(5) = 0.236926885056189087514264040720;
        
    otherwise        
        error ( 'Invalid number of quad points.' );
end




