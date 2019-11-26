% *************************************************************************
% Filename: MAE-5700_HW11.m
% Author: Dalton and Sairam
% Date: 11/22/19
% Notes: Reduced code modifcations for HW11 2DBVP for heat convection
% *************************************************************************

% From InputData.m ***************************************************
% ...

% point source boundary condition for heat convection, note that there must
% be a node at the specific location of the point source application
% [x-coord, y-coord, heat source value]
boundStruct.PointSource = [0 5 125];


% From ApplyBC.m ***************************************************
% ...

% Apply point source, collect necessary variables from structs
PointSource = boundStruct.PointSource;
nCoords = meshStruct.nCoords;

% Iterate over each point source 
for i = 1:size( boundStruct.PointSource,1
   % get the current point source
   curr_PS =  boundStruct.PointSource(i,:);
   
   % find the node number associated with the coordinates of the point
   % source
   nodeNum_PS = find(meshStruct.nCoords(:,1) == curr_PS(1) & meshStruct.nCoords(:,2) == curr_PS(2));
   
   % apply point source to force vector at node location
   globalSystem.F(nodeNum_PS) = globalSystem.F(nodeNum_PS) + curr_PS(3);
end
