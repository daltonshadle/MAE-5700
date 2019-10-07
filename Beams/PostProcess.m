function globalSystem=PostProcess(globalSystem,meshStruct,boundStruct);
% globalSystem=POSTPROCESS(globalSystem,meshStruct,boundStruct)
% Calculate bending moment and shear force for BEAM problem.
% Plot the deformed shapes of the beam.
% Print relevant problem data and results.
% last edit: 5 August 2015 H. Ritz

% unpack necessary input
numEls    =meshStruct.numEls;
gatherMat =meshStruct.gatherMat;
elEI      =meshStruct.elEI;
nCoords   =meshStruct.nCoords;
elCon     =meshStruct.elCon;
d         =globalSystem.d;

% initialize vectors for plotting
displacement = [];
moment       = [];
shear        = [];

for elmID = 1:numEls % loop over elements to plot
                     % displacements, moments and shear forces
                     
    elmSoln=d(gatherMat(elmID,:)); % extract element nodal displacements
    EI   = elEI(elmID);            % Young's modulus x I
    
    xe = nCoords(elCon(elmID,:));   % extract element coordinates
    L  = xe(2) - xe(1); % element length
    J  = L / 2;         % compute Jacobian for transformation.
    
    nplot = 10;                  % number of points to plot
    
    % Compute displacements, moments and shear forces
    x=linspace(xe(1),xe(2),nplot); % equally distributed points in element
    xplot = 2/L*(x-xe(1))-1;       % transform to local coordinate
    
    for i = 1:nplot
        xi = xplot(i);    % current coordinate
        % shape functions evaluated at this point
        N = [1/4*(1-xi)^2*(2+xi), L/8*(1-xi)^2*(1+xi), ...
            1/4*(1+xi)^2*(2-xi), L/8*(1+xi)^2*(xi-1)];  
        % 2nd derivative of N
        B = [3/2*xi, L*(3/4*xi - 1/4), -3/2*xi,L*(3/4*xi + 1/4)]/J^2;
        % 3rd derivative of N
        S = [3/2, 3/4*L, -3/2, 3/4*L]/J^3;
        
        % The following commands append as a new row the vector consisting
        % of x location & corresponding displacement/moment/shear force.
        % This allows easy plotting of the displacements, moments and shear
        % force as a function of location in the beam.
        
        displacement    = [displacement;[x(i),  N*elmSoln]]; 
        moment          = [moment; [x(i), EI*B*elmSoln]];    
        shear           = [shear;  [x(i), -EI*S*elmSoln]];   
    end
end
% % save a workspace with certain outputs
% eval(['displacement',num2str(numEls),'=displacement;'])
% eval(['moment',num2str(numEls),'=moment;'])
% eval(['shear',num2str(numEls),'=shear;'])
% outputfile=sprintf('Output%d.mat',numEls)
% eval(['save(outputfile, ''displacement',num2str(numEls),''', ''shear',num2str(numEls),''', ''moment',num2str(numEls),''');'])

% break
% plot displacements, moment and shear forces
set(0,'defaultLineLineWidth',3)
set(0,'defaultTextFontSize',24)
set(0,'defaultAxesFontSize',24)
set(0,'defaultAxesFontWeight','bold')

% calculate exact solution
[x_plot,S_ex,M_ex,w_ex]=Exact;
figure
plot(displacement(:,1),displacement(:,2),'-r'); hold on;
plot(x_plot,w_ex,'-k'); legend('FE','Exact Solution'); hold on;
ylabel('displacement');  xlabel('x'); title(['Displacements (',num2str(numEls),' elements)']);

% break
figure
plot(moment(:,1),moment(:,2),'-r'); hold on;
plot(x_plot,M_ex,'-k'); legend('FE','Exact Solution'); hold on;
ylabel('moment'); xlabel('x'); title(['Moments (',num2str(numEls),' elements)']);

figure
plot(shear(:,1),shear(:,2),'-r'); hold on;
plot(x_plot,S_ex,'-k'); legend('FE','Exact Solution'); hold on;
ylabel('shear'); xlabel('x');   title(['Shear (',num2str(numEls),' elements)']);

% Package variables into the output struct
globalSystem.displacement=displacement;
globalSystem.moment      =moment;
globalSystem.shear       =shear;

