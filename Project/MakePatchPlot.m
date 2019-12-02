function MakePatchPlot(z,range,meshstruct,titlestring)
% MakePatchPlot(Z,RANGE,MESHSTRUCT) 
%
% This function accepts a vector Z, a scalar RANGE, and a structure MESHSTRUCT including
% coordinates and connectivity. It then plots Z in the current axes on 
% an element-by-element basis (patch plot). RANGE is used to determine the
% limits on the colorbar.

% initialize variables
% % meshdata=meshstruct.size;
numel = meshstruct.numEls;
nnpe = meshstruct.nnpe;
coords=meshstruct.nCoords;
x=coords(:,1);
y=coords(:,2);
np=meshstruct.elCon;

% find appropriate colorbar range
cmax = max(z);
cmin = min(z);
% if abs(cmax-cmin)<(range*.001)
%     tmp=0.10*range;
%     cmax=cmax+tmp
%     cmin=cmin-tmp
% end

% prepare axes
figure;
a = [cmin cmax];
caxis(a);
axis equal;
hold on

% plot one patch at a time
for  i=1:numel
    for j=1:nnpe
        j1=np(i,j);
        xvertex(j) = x(j1);
        yvertex(j) = y(j1);
        zval(j)    = z(j1);
    end
    cval=zval;
%     cval=(zval-cmin)./(cmax-cmin)
    patch(xvertex',yvertex',cval','EdgeColor','none','FaceColor','interp');
end
hold off
title(titlestring);
h=colorbar('EastOutside');
if ~ispc
    set(h,'FontSize',14);
end
