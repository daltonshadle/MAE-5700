function errorStruct=ErrorAnalysis(globalSystem,meshStruct,ufighand,dufighand)
% errorStruct=ERRORANALYSIS(globalSystem,meshStruct,uFEM,duFEM,ufighand,dufighand) 
% This compares the finite element solution and derivative to the exact
% solution and derivative. It also calculates the L_2 error nrom and the
% energy error norm. Only call this script if the exact solution is
% available for the current problem.
% last edit: 20 October 2015 H. Ritz

% unpack necessary inputs
gatherMat=meshStruct.gatherMat;
numEls   =meshStruct.numEls;
numQP    =meshStruct.numQP;
nCoords  =meshStruct.nCoords;
nnpe     =meshStruct.nnpe;
d        =globalSystem.d;

% initialize variables
energyNorm = 0;  % Energy norm, from derivatives at quad points
L2Norm = 0.0;    % L2 norm, from solution at quad points

energyNormDenom = 0; % Normalization factor for the energy
                     % norm error (based on the exact
                     % solution)

L2NormDenom = 0;     % Normalization factor for the
                     % L2 error norm (based on the exact
                     % solution)

% extract Gauss points and weights
[qp,w] = Gauss(5); % use the most accurate integration possible

for elmID = 1:numEls % loop over elements
    gdof = gatherMat(elmID,:); % Extract the global degrees of freedom
    xvec = nCoords(gdof)';     % Extract the coordinate vector for this element
    uvec = d(gdof);            % Extract the solution vector for this element
    
    for iqp = 1:length(qp)   % loop over all the quad points
        NofXi=Nmatrix(qp(iqp),nnpe); % row vector of shape functions
                                % at this XI location in parent domain
        dNdXi=Bmatrix(qp(iqp),nnpe); % row vector of shape function derivatives
                                % at this XI location in parent domain
        JofXi=dNdXi*xvec;  % jacobian at this XI location (scalar)
        B=dNdXi/JofXi;     % row vector of dNdX at this XI location
        
        X=NofXi*xvec;      % physical coordinate of this XI location (scalar)
        
        [u du] = Exact(X); % Analytical solution
        
        UatXi  = NofXi*uvec; % solution at XI
        dUatXi = B*uvec;     % derivative at XI
        
        uerror =abs(UatXi - u); % error at each quad point
        duerror=abs(dUatXi - du); % error at each quad point
        L2Norm     = L2Norm + (uerror)^2*JofXi*w(iqp);
        L2NormDenom =  L2NormDenom + u^2*JofXi*w(iqp);
        
        energyNorm = energyNorm + 1/2*PP(X)*(duerror)^2*JofXi*w(iqp);
        energyNormDenom = energyNormDenom + 1/2*PP(X)*du^2*JofXi*w(iqp);
    end
end

L2Norm = sqrt(L2Norm) / sqrt (L2NormDenom);

energyNorm = sqrt(energyNorm) / sqrt (energyNormDenom) ;

% package outputs
errorStruct.L2Norm=L2Norm;
errorStruct.energyNorm=energyNorm;

fprintf(1,'\nError Reports (%d elements):\n',numEls);
fprintf(1,'\tThe normalized L_2 norm error is %e\n', L2Norm);
fprintf(1,'\tThe normalized energy norm error is %e\n', energyNorm);

x_plot = linspace (nCoords(1), nCoords(end), 1000); % Extract analytical
[u du] = Exact(x_plot);                        % solution

% Plot the solution
figure(ufighand)
plot(x_plot,u,'-k'); legend('FE','Exact Solution','Location','Best'); hold on;
title(['u: FE (',num2str(numEls),' els) vs. analytical']);

figure(dufighand)
plot(x_plot,du,'-k'); legend('FE','Exact Solution','Location','Best'); hold on;
title(['du/dx: FE (',num2str(numEls),' els) vs. analytical']);






