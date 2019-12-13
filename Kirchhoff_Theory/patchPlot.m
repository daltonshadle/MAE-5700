function patchPlot(fieldValue, meshStruct, plotTitle, filename, varargin)
%patchPlot(fieldValue, meshStruct, plotTitle,[nColor])        contour plot
% 
% Input
% fieldValue : column or row vector of discrete field values for contour
% plotting;
% meshStruct : structure for mesh info. See LinElast.m for detailed specs.
% plotTitle : string for the figure title
% nColor (optional arguement) : how many contours to use

% Last editted Nov. 22, 2017, H. Ritz

% number of different colors in the contour plot
if isempty(varargin)
    nColor  = 10;
else
    nColor=varargin{1};
end

% make fieldValue a column vector
if isrow(fieldValue)
    fieldValue = fieldValue';
end

nCoords=meshStruct.nCoords;
elCon=meshStruct.elCon;

% xTri (yTri) is the x (y) coordinates of the nodals of the triangulars
xCoords = nCoords(:,1);
yCoords = nCoords(:,2);
xElm = xCoords(elCon)'; 
yElm = yCoords(elCon)';

% field values at the nodals of the triangulars
fieldElm = fieldValue(elCon)';

% plotting
gcf = figure;
x0=10;
y0=10;
width=800;
height=700;
set(gcf,'position',[x0,y0,width,height])
patch(xElm,yElm,fieldElm,'EdgeColor','interp');

colormap(jet(nColor));
hColorbar = colorbar;

% Tick markers positioned in between color transitions
if nColor<=15
    nTick = nColor+1;
else
    nTick = 16;
end
hColorbar.Ticks = linspace(min(fieldValue),max(fieldValue),nTick);
xlabel('x');
ylabel('y');

title(plotTitle);
axis image

saveas(gcf,filename,'epsc')
saveas(gcf,filename,'png')