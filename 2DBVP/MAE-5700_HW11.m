% *************************************************************************
% Filename: MAE-5700_HW11.m
% Author: Dalton and Sairam
% Date: 11/22/19
% Notes: Reduced code modifcations for HW11 2DBVP for torsion
% *************************************************************************

% From PostProccessor.m ***************************************************
% ...

% Pull necessary variables from structs
d = globalSystem.d;
nnpe = meshStruct.nnpe;
numEls = meshStruct.numEls;
elCon = meshStruct.elCon;
nCoords = meshStruct.nCoords;

% Calculate gradient
gradient = calNodalGradient(d,meshStruct);

% Get stresses
sigma_13 = gradient(:,2);
sigma_23 = -gradient(:,1);

% Plot Stress Vector
if strcmp(PlotInstructions.plot_vector,'yes')
    plotvector([sigma_13,sigma_23],meshStruct.nCoords,'FE solution: Stress Vectors',1);
end

% Plot the contour distribution of stress fields
if strcmp(PlotInstructions.plot_contour,'yes')
    patchPlot(sigma_13, meshStruct,'\sigma_{13} field',15);
end
% Plot the contour distribution of stress fields
if strcmp(PlotInstructions.plot_contour,'yes')
    patchPlot(sigma_23, meshStruct,'\sigma_{23} field',15);
end


% Calc von Mises effective stress
von_Mises = sqrt(sigma_13.^2 + sigma_23.^2);

% Plot the contour of von Mises
if strcmp(PlotInstructions.plot_contour,'yes')
    patchPlot(von_Mises, meshStruct,'von Mises field',15);
end

% Calculate torque
tot_torque = 0;
switch nnpe
    case 3 % Triangle Elements
        % iterate over each element
        for i = 1:numEls
            % collect nodal coordinates
            elNodes = elCon(i,:);
            elCoords = nCoords(elNodes,:);
            a = elCoords(1,:);
            b = elCoords(2,:);
            c = elCoords(3,:);
            
            % calculate area with coordinates
            elArea = (a(1)*(b(2)-c(2)) + b(1)*(c(2)-a(2)) + c(1)*(a(2)-b(2)))/2;
            
            % average phi from nodes (Note: we're assuming we can simply
            % average the nodal values since our elements are small enough
            % that error won't be affected to drastically)
            elPhi = sum(d(elNodes)) / nnpe;
            
            % add elemental torque to total
            tot_torque = tot_torque + 2 * elArea * elPhi;
        end
    case 4 % Quad Elements
        % iterate over each element
        for i = 1:numEls
            % collect nodal coordinates
            elNodes = elCon(i,:);
            elCoords = nCoords(elNodes,:);
            a = elCoords(1,:);
            b = elCoords(2,:);
            c = elCoords(3,:);
            e = elCoords(4,:);
            
            % calculate area with coordinates
            elArea = ((a(1)*b(2) - a(2)*b(1)) + (b(1)*c(2) - b(2)*c(1)) ...
                + (c(1)*e(2) - c(2)*e(1)) + + (e(1)*a(2) - e(2)*a(1)))/2;
            
            % average phi from nodes (Note: we're assuming we can simply
            % average the nodal values since our elements are small enough
            % that error won't be affected to drastically)
            elPhi = sum(d(elNodes)) / nnpe;
            
            % add elemental torque to total
            tot_torque = tot_torque + 2 * elArea * elPhi;
        end
end

% Display total torque
fprintf('Total Torque: %.3f\n', tot_torque);
