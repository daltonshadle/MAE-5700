%% Jacobian variation in Quadrilateral elements
% Use this function to explore the value of the Jacobian at different
% locations in an element based on nodal locations.
function HW7Problem2

%% Set up

clear 
close all
set(0,'defaultLineLineWidth',3)
set(0,'defaultTextFontSize',24)
set(0,'defaultAxesFontSize',24)
set(0,'defaultAxesFontWeight','bold')

%% Q4 elements
% For the assigned nodal positions, show that the Jacobian can be negative
% for concave Q4 elements.

%%
% Input global nodal positions
XY=[2 1 ; 0 0; 1 0; 2 -1];
JacCalcFun(XY)
%% Q8 elements
% See how the Jacobian varies as the midside node moves for a Q8 element.
%%
% Input global nodal positions
XY=[0 0; 1 0 ; 1 1 ; 0 1 ; 0.5 0 ; 1 0.5 ; 0.5 1 ; 0 0.5];
JacCalcFun(XY);

%% Calculation function
% This function actually contains all calculations, based on the nodal
% positions.
function JacCalcFun(XY)
nnpe=size(XY,1);
%%
% Plot the element shape. For 8 node elements, plot the curved sides.
figure
switch nnpe
    case 4
        plot(XY([1:end,1],1),XY([1:end,1],2),'o-'); axis equal
    case 8
        edges=linspace(-1,1,10);
        exi=[edges,ones(size(edges)),-edges,-1*ones(size(edges))];
        eeta=[-1*ones(size(edges)),edges,ones(size(edges)),-edges];
        for e=1:length(exi)
            N=Nmat(exi(e),eeta(e),nnpe); % get the shape functions to plot in the physical domain
            ex(e)=N*XY(:,1); % use the shape functions to find the global
            ey(e)=N*XY(:,2); % coordinates of this \xi, \eta
        end
        plot(ex,ey,'o-'); axis equal
end
        


%% 
% Make a grid of $\xi$ and $\eta$ positions to calculate the Jacobian. Use
% the shape functions to find the physical coordinate of the $(\xi, \eta)$
% position.
points=linspace(-1,1,100);
[XI,ETA]=meshgrid(points,points);

for row=1:size(XI,1)
    for col=1:size(XI,2)
        xi=XI(row,col); eta=ETA(row,col);
        Jac=JacMat(xi,eta,XY);
        detJ(row,col)=det(Jac);
        N=Nmat(xi,eta,nnpe); % get the shape functions to plot in the physical domain
        x(row,col)=error('Find the x position for this point using N.') % use the shape functions to find the global 
        y(row,col)=error('Find the y position for this point using N.') % coordinates of this \xi, \eta
    end
end

%%
% Plot the determinant of the Jacobian in both parent and physical domains.

figure
subplot(1,2,1)
surf(XI,ETA,detJ,'EdgeColor','none')
view(2)
colorbar
xlabel('\xi')
ylabel('\eta')
axis square

subplot(1,2,2)
surf(x,y,detJ,'EdgeColor','none')
view(2)
colorbar
xlabel('x')
ylabel('y')
axis equal tight



%% Shape function array
function N=Nmat(xi, eta, nnpe)
switch nnpe
    case 4 
        N=1/4*[(1-xi)*(1-eta), (1+xi)*(1-eta), (1+xi)*(1+eta), (1-xi)*(1+eta)];
    case 8
        N=error('Define the shape function array for Q8 elements.')
end        

%% Shape function gradient array
function GN=GNmat(xi, eta, nnpe)
switch nnpe
    case 4
        GN=error('Define the shape function gradient array for Q4 elements.')
    case 8
        GN=[
            -1/4*(-1+eta)*(+2*xi+eta),-1/4*(-1+xi)*(+xi+2*eta), 
            +1/4*(-1+eta)*(-2*xi+eta),+1/4*(+1+xi)*(-xi+2*eta), 
            +1/4*(+1+eta)*(+2*xi+eta),+1/4*(+1+xi)*(+xi+2*eta), 
            -1/4*(+1+eta)*(-2*xi+eta),-1/4*(-1+xi)*(-xi+2*eta), 
            xi*(-1+eta),1/2*(1+xi)*(-1+xi),
            -1/2*(1+eta)*(-1+eta),-eta*(1+xi),
            -xi*(1+eta),-1/2*(1+xi)*(-1+xi),
            1/2*(1+eta)*(-1+eta),eta*(-1+xi),
            ]'
end        

%% Jacobian matrix
function Jac=JacMat(xi, eta, XY);
GN=GNmat(xi,eta,size(XY,1));
Jac=error('Calculate the Jacobian matrix at this location.')
