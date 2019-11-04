function globalSystem=PostProcess(globalSystem,meshStruct)
% This script plots the solution and the derivative of the solution, and
% analyzes the error if an exact solution is available.
% last edit: 6 August 2015 H. Ritz

% unpack necessary inputs
gatherMat=meshStruct.gatherMat;
numEls   =meshStruct.numEls;
nCoords  =meshStruct.nCoords;
nnpe     =meshStruct.nnpe;
d        =globalSystem.d;

% initialize the vectors for plotting
u_plot =[]; du_plot =[];

nplot  = 10;    % number of plot points per element
xiplot = linspace(-1,1,nplot); % xi values to plot

for elmID = 1:numEls
    gdof = gatherMat(elmID,:); % Extract the global degrees of freedom
    xvec = nCoords(gdof)';     % Extract the coordinate vector for this element
    uvec = d(gdof);            % Extract the solution vector for this element
    
    for ixi = 1:nplot          % loop over all the plotting points
        xi=xiplot(ixi);    % local coordinate
        NofXi=Nmatrix(xi,nnpe); % row vector of shape functions
                           % at this XI location in parent domain
        dNdXi=Bmatrix(xi,nnpe); % row vector of shape function derivatives
                           % at this XI location in parent domain
        JofXi=dNdXi*xvec;  % jacobian at this XI location (scalar)
        B=dNdXi/JofXi;     % row vector of dNdX at this XI location
        
        X=NofXi*xvec;      % physical coordinate of this XI location (scalar)

        UatXi  = NofXi*uvec; % solution at XI
        dUatXi = B*uvec;     % derivative at XI
        
        u_plot  = [u_plot; X UatXi];
        du_plot = [du_plot; X dUatXi];
    end
end
% Save output
globalSystem.uPlot=u_plot;
globalSystem.duPlot=du_plot;

% set defaults for plots
set(0,'defaultLineLineWidth',3)
set(0,'defaultTextFontSize',24)
set(0,'defaultAxesFontSize',24)
set(0,'defaultAxesFontWeight','bold')

fighandsoln=figure;  % Plot the FE solution
plot(u_plot(:,1),u_plot(:,2),'r');
hold on;
ylabel('u');  xlabel('x'); title(['u: FE solution (',num2str(numEls),' els)']);

fighandderiv=figure;  % Plot the derivative of the solution
plot(du_plot(:,1),du_plot(:,2),'r');
hold on;
ylabel('du/dx'); xlabel('x'); title(['du/dx: FE solution (',num2str(numEls),' els)']);



% If you have an exact solution for this problem, use the following code.
% If not, comment this section out of the code.
globalSystem.errorStruct=ErrorAnalysis(globalSystem,...
    meshStruct,fighandsoln,fighandderiv); % plot the exact solution and 
                                          % calculate the error
