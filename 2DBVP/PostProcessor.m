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

if_torsion = 0;

% Plot the contour distribution
if strcmp(PlotInstructions.plot_contour,'yes')
    patchPlot(d, meshStruct,'Temperature Field',15);
end

% Calculate gradient
gradient = calNodalGradient(d,meshStruct);

if strcmp(PlotInstructions.plot_vector,'yes')
    plotvector(-gradient*DD(1),meshStruct.nCoords,'FE solution: Flux',1);
end


if (if_torsion)
    % Get stresses
    sigma_13 = gradient(:,2);
    sigma_23 = -gradient(:,1);
    
    % Plot Traction Vector
    if strcmp(PlotInstructions.plot_vector,'yes')
        plotvector([sigma_13,sigma_23],meshStruct.nCoords,'FE solution: Stress Vectors',1);
    end
    % Plot the contour distribution
    if strcmp(PlotInstructions.plot_contour,'yes')
        patchPlot(sigma_13, meshStruct,'\sigma_{13} field',15);
    end% Plot the contour distribution
    if strcmp(PlotInstructions.plot_contour,'yes')
        patchPlot(sigma_23, meshStruct,'\sigma_{23} field',15);
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
    
    
    fprintf('Total Torque: %.3f\n', tot_torque);
end
