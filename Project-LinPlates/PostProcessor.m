% *************************************************************************
% Name: PostPRocessor.m
% Notes: PostProcessor(PlotInstructions,meshStruct,globalSystem)
%     Calculate stress, make contour plots, vector field plots, deformed 
%     mesh, etc. 
%     last edit: 22 November 2017 H. Ritz
% Project Updates: Pulled from LinElast code base from MAE-5700 course,
%     updates were made to implement Kirchhoff Theory Plate Bending
% Update Authors: Dalton and Sairam
% *************************************************************************

function [globalSystem] = PostProcessor(PlotInstructions,meshStruct,globalSystem)
%% Unpack variables from structures
d = globalSystem.d;
neq=meshStruct.numEq;
nCoords = meshStruct.nCoords;
elCon = meshStruct.elCon;
thickness_layer = meshStruct.thickness_layer;

%% capturing elemental displacements and angular strains
u_z = d(1:3:neq-2);
theta_x = d(2:3:neq-1);
theta_y = d(3:3:neq);

%% gca properties
set(0,'defaultLineLineWidth',1.5)
set(0,'defaultTextFontSize',12)
set(0,'defaultAxesFontSize',12)
set(0,'defaultAxesFontWeight','bold')

%% Plot Deformed Mesh
loadStr = 'vary_clamped_';
meshStr = [num2str(meshStruct.nx), 'x' num2str(meshStruct.ny)];
if strcmp(PlotInstructions.plot_deformed,'yes')
    PlotDeformedMesh('Deformed Configuration',meshStruct,d, [loadStr, 'deform', meshStr]);
end

%% Stresses and strains
top  = length(thickness_layer);
thickness_layer(top)
bottom = 1;
thickness_layer(bottom)
mid = abs(bottom - top)/2 +1;

[strain,stress,moment] = GPNodalStrainStress(d,meshStruct,thickness_layer(bottom));
strain = full(strain); stress = full(stress); moment = full(moment);
globalSystem.strain = strain;
globalSystem.stress = stress;
globalSystem.moment = moment;

% plotting stresses
patchPlot(stress(:,1), meshStruct, '\sigma_{xx} (Pa)', [loadStr, 'sigma_xx_', meshStr]);
patchPlot(stress(:,2), meshStruct, '\sigma_{yy} (Pa)', [loadStr, 'sigma_yy_', meshStr]);
patchPlot(stress(:,3), meshStruct, '\sigma_{xy} (Pa)', [loadStr, 'sigma_xy_', meshStr]);

vM_stress = sqrt(1/2 * ((stress(:,1) - stress(:,2)).^2 + stress(:,1).^2 + stress(:,2).^2 + 6*stress(:,3).^2));
patchPlot(vM_stress, meshStruct, 'von-Mises Stress (\sigma_v) (Pa)', [loadStr, 'vM_stress_', meshStr]);

% plotting strains
patchPlot(strain(:,1), meshStruct, '\epsilon_{xx}', [loadStr, 'epsilon_xx_', meshStr]);
patchPlot(strain(:,2), meshStruct, '\epsilon_{yy}', [loadStr, 'epsilon_yy_', meshStr]);
patchPlot(strain(:,3), meshStruct, '\epsilon_{xy}', [loadStr, 'epsilon_xy_', meshStr]);

% plotting moments
patchPlot(moment(:,1), meshStruct, 'Moment per unit length, M_x (N)', [loadStr, 'M_x_', meshStr]);
patchPlot(moment(:,2), meshStruct, 'Moment per unit length, M_y (N)', [loadStr, 'M_y_', meshStr]);
patchPlot(moment(:,3), meshStruct, 'Moment per unit length, M_{xy} (N)', [loadStr, 'M_xy_', meshStr]);
















%%

% if strcmp(PlotInstructions.plot_deformed,'yes')
%     PlotDeformedMesh('Deformed Configuration',meshStruct,d);
% end

% %% Plot the contour distribution
% if strcmp(PlotInstructions.plot_contour,'yes')
%     patchPlot(u,meshStruct, 'FE solution: u displacement ');
%     patchPlot(v,meshStruct, 'FE solution: v displacement ');
% end
%
% Calculate the strain and stress


% sigma_xx = sig(:,1);
% sigma_yy = sig(:,2);
% sigma_xy = sig(:,3);
%
% % Calculate the principal stress
% sigma_1 = (sigma_xx + sigma_yy)/2 + sqrt(((sigma_xx + sigma_yy)/2).^2 + sigma_xy.^2);
% sigma_2 = (sigma_xx + sigma_yy)/2 - sqrt(((sigma_xx - sigma_yy)/2).^2 + sigma_xy.^2);
% tau_max = (sigma_1-sigma_2)/2;
% sigma_vm=sqrt(sigma_xx.^2+sigma_yy.^2-sigma_xx.*sigma_yy+3*sigma_xy.^2);
%
% % Plot stresses
% if strcmp(PlotInstructions.plot_contour,'yes')
%
%     patchPlot(sigma_xx,meshStruct, 'xx stress ');
%     patchPlot(sigma_yy,meshStruct, 'yy stress ');
%     patchPlot(sigma_xy,meshStruct, 'xy stress ');
%     patchPlot(sigma_vm,meshStruct, 'von Mises stress ');
% end
%
% if strcmp(PlotInstructions.plot_fringes,'yes')
%     PlotFringes(2*tau_max, 10,meshStruct);
% end
