function PlotFringes(tau,numfringe,MeshStruct)
% PlotFringes
% This function  plots the fringe pattern
% Input tau -- maximum in plane shear stress (sigma_1 - sigma_2)
% Input numfringe -- number of fringes

Nodes=MeshStruct.nCoords';
Elems=MeshStruct.elCon';
nno=MeshStruct.numNodes;
nen=MeshStruct.nnpe;
nel=MeshStruct.numEls;

figure;
% make fringes
high_shear = max(tau);
low_shear = min(tau);
shear_span = high_shear-low_shear;

if(high_shear == 0.0)
    error(' No fringes ');
end

if(shear_span < 0.02*high_shear)
    numfringe = 1.0;
    shear_span = high_shear;
end

fringe_vect = (numfringe/shear_span)*(tau - low_shear*ones(nno,1));
dark_values = fix(fringe_vect);
fringe_shade = sin((fringe_vect - dark_values)*pi);

% use grayscale colormap
colormap(gray);
cmax = 1.0; cmin = 0.0; a = [cmin cmax]; caxis(a);
for i=1:nel

    glb = Elems(:,i)';           %Extract the global node numbers

    XX = [Nodes(1, glb),Nodes(1,glb(1))];
    YY = [Nodes(2, glb),Nodes(2,glb(1))];

    dd = [fringe_shade(glb)', fringe_shade(glb(1))];

    patch(XX,YY,dd,'EdgeColor','interp','FaceColor','interp');

    hold on;

end
axis image;
title('Fringe pattern');

