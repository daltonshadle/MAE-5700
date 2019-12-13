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
 function PlotDeformedMesh(string_title,meshStruct,d)
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

nel = size(elCon,1) ;              % number of elements
nnode = size(nCoords,1) ;          % total number of nodes in system
nnpe = size(elCon,2);              % number of nodes per element

% Initialization of the required matrices
X = zeros(nnpe,nel) ;
Y = zeros(nnpe,nel) ;
Z = zeros(nnpe,nel) ;
profile = zeros(nnpe,nel) ;

% top layer % undeformed
for iel=1:nel % loop over elements
    for i=1:nnpe % number of nodes per element
        nd(i)=elCon(iel,i);         % extract connected node for (iel)-th element
        X(i,iel)=nCoords(nd(i),1);    % extract x value of the node
        Y(i,iel)=nCoords(nd(i),2);    % extract y value of the node
    end
    Z(:,iel) = thickness_layer(end)*ones(4,1);
    profile(:,iel) = 0;
end

% top layer % deformed
for iel=1:nel % loop over elements
    for i=1:nnpe % number of nodes per element
        nd(i)=elCon(iel,i);         % extract connected node for (iel)-th element
        X_def(i,iel)= -thickness_layer(end)*theta_x(nd(i))+ nCoords(nd(i),1);    % extract x value of the node
        Y_def(i,iel)= -thickness_layer(end)*theta_x(nd(i))+ nCoords(nd(i),2);    % extract y value of the node
    end
    Z_def(:,iel) = u_z(nd)' + thickness_layer(end)*ones(4,1);
    profile_def(:,iel) = u_z(nd');
end


% Plotting the profile of a property on the deformed mesh
fh = figure ;
subplot(2,1,1)
fill3(X,Y,Z,profile)
zlim([thickness_layer(1) thickness_layer(end)])
title('Undeformed configuration')
subplot(2,1,2)
set(fh,'name','Postprocessing','numbertitle','off') ;
fill3(X_def,Y_def,Z_def,profile_def)
title(string_title)
% Colorbar Setting
SetColorbar

