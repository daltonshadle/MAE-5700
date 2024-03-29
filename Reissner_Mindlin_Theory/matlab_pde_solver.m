clc;
clear all;
close all;
%%
numberOfPDE = 2;
model = createpde(numberOfPDE);
%%
E = 1.0e6; % modulus of elasticity
nu = .3; % Poisson's ratio
thick = .1; % plate thickness
len = 10.0; % side length for the square plate
hmax = len/20; % mesh size parameter
D = E*thick^3/(12*(1 - nu^2));
pres = 2; % external pressure
%%
gdm = [3 4 0 len len 0 0 0 len len]';
g = decsg(gdm,'S1',('S1')');
%%
geometryFromEdges(model,g);
%%
figure; 
pdegplot(model,'EdgeLabels','on');
ylim([-1,11])
axis equal
title 'Geometry With Edge Labels Displayed';
%%
c = [1 0 1 D 0 D]';
a = [0 0 1 0]';
f = [0 pres]';
specifyCoefficients(model,'m',0,'d',0,'c',c,'a',a,'f',f);
%%
k = 1e7; % spring stiffness
%%
bOuter = applyBoundaryCondition(model,'neumann','Edge',(1:4),...
                                     'g',[0 0],'q',[0 0; k 0]);
%%
generateMesh(model, 'Hmax', hmax);
%%
res = solvepde(model);
u = res.NodalSolution;
numNodes = size(model.Mesh.Nodes,2);
figure
pdeplot(model,'XYData',u(1:numNodes),'Contour','on');
title 'Transverse Deflection'
%%

numNodes = size(model.Mesh.Nodes,2);
fprintf('Transverse deflection at plate center(PDE Toolbox) = %12.4e\n',...
                                                  min(u(1:numNodes,1)));
 %%
wMax = -.0138*pres*len^4/(E*thick^3);
fprintf('Transverse deflection at plate center(analytical) = %12.4e\n', wMax);
%%





