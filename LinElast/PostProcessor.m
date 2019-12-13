% PostProcessor(PlotInstructions,meshStruct,globalSystem)
% Calculate stress, make contour plots, vector field plots, deformed mesh, etc.
%
% last update: 22 November 2017 H. Ritz
function [globalSystem] = PostProcessor(PlotInstructions,meshStruct,globalSystem)
d = globalSystem.d;
neq=meshStruct.numEq;
u = d(1:2:neq-1);
v = d(2:2:neq);

set(0,'defaultLineLineWidth',1.5)
set(0,'defaultTextFontSize',12)
set(0,'defaultAxesFontSize',12)
% set(0,'defaultAxesFontWeight','bold')


if strcmp(PlotInstructions.plot_deformed,'yes')
    PlotDeformedMesh('Deformed Configuration',meshStruct,d);
end
% Plot the contour distribution
% if strcmp(PlotInstructions.plot_contour,'yes')
%     patchPlot(u,meshStruct, 'FE solution: u displacement ');
%     patchPlot(v,meshStruct, 'FE solution: v displacement ');
% end

% Calculate the strain and stress
[~,sig] = calNodalStrainStress(d,meshStruct);

sigma_xx = sig(:,1);
sigma_yy = sig(:,2); 
sigma_xy = sig(:,3);

% Calculate the principal stress
sigma_1 = (sigma_xx + sigma_yy)/2 + sqrt(((sigma_xx + sigma_yy)/2).^2 + sigma_xy.^2);
sigma_2 = (sigma_xx + sigma_yy)/2 - sqrt(((sigma_xx - sigma_yy)/2).^2 + sigma_xy.^2);
tau_max = (sigma_1-sigma_2)/2;
sigma_vm=sqrt(sigma_xx.^2+sigma_yy.^2-sigma_xx.*sigma_yy+3*sigma_xy.^2);

% Plot stresses
if strcmp(PlotInstructions.plot_contour,'yes')

    patchPlot(sigma_xx,meshStruct, '\sigma_{xx}');
    patchPlot(sigma_yy,meshStruct, '\sigma_{yy}');
    patchPlot(sigma_xy,meshStruct, '\sigma_{xy}');
    patchPlot(sigma_vm,meshStruct, 'von Mises stress ');
end

if strcmp(PlotInstructions.plot_fringes,'yes')
    PlotFringes(2*tau_max, 10,meshStruct);
end

globalSystem.sig = sig;
globalSystem.v = v;