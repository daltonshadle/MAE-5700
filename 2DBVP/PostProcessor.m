function PostProcessor(PlotInstructions,meshStruct,globalSystem)
% PostProcessor(PlotInstructions,meshStruct,globalSystem)
% Calculate flux, make contour plots, vector field plots, etc.
%
% last update: 21 November 2017 H. Ritz

d = globalSystem.d;
nnpe = meshStruct.nnpe;
numEls = meshStruct.numEls;
elCon = meshStruct.elCon;
nCoords = meshStruct.nCoords;

% Plot the contour distribution
if strcmp(PlotInstructions.plot_contour,'yes')
    patchPlot(d, meshStruct,'\phi field',15);
end

% Calculate gradient
gradient = calNodalGradient(d,meshStruct);

% Plot Gradient
if strcmp(PlotInstructions.plot_vector,'yes')
    plotvector(gradient,meshStruct.nCoords,'FE solution: \nabla\phi',1);
end

% Get stresses
sigma_13 = gradient(:,2);
sigma_23 = -gradient(:,1);

% Plot Traction Vector
if strcmp(PlotInstructions.plot_vector,'yes')
    plotvector([sigma_13,sigma_23],meshStruct.nCoords,'FE solution: Traction Vectors',1);
end

% Calc von Mises effective stress
von_Mises = sqrt(sigma_13.^2 + sigma_23.^2);

% Plot the contour of von Mises
if strcmp(PlotInstructions.plot_contour,'yes')
    patchPlot(von_Mises, meshStruct,'von Mises field',15);
end


switch nnpe
    case 3 % Triangle Elements
        tot_torque = 0;
        for i = 1:numEls
           elNodes = elCon(i,:);
           elCoords = nCoords(elNodes,:);
           a = elCoords(1,:);
           b = elCoords(2,:);
           c = elCoords(3,:);

           elArea = (a(1)*(b(2)-c(2)) + b(1)*(c(2)-a(2)) + c(1)*(a(2)-b(2)))/2;

           elPhi = sum(d(elNodes)) / nnpe;

           tot_torque = tot_torque + 2 * elArea * elPhi;
        end
    case 4 % Quad Elements
        tot_torque = 0;
        for i = 1:numEls
           elNodes = elCon(i,:);
           elCoords = nCoords(elNodes,:);
           a = elCoords(1,:);
           b = elCoords(2,:);
           c = elCoords(3,:);
           e = elCoords(4,:);
           
           elArea = ((a(1)*b(2) - a(2)*b(1)) + (b(1)*c(2) - b(2)*c(1)) ...
               + (c(1)*e(2) - c(2)*e(1)) + + (e(1)*a(2) - e(2)*a(1)))/2;
           
           elPhi = sum(d(elNodes)) / nnpe;
           
           tot_torque = tot_torque + 2 * elArea * elPhi;
        end
end


disp(tot_torque);
