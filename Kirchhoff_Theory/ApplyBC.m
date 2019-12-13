% *************************************************************************
% Name: ApplyBC.m
% Notes: [globalSystem,boundStruct] = ApplyBC(boundStruct,meshStruct,globalSystem)
%     This function applies natural and essential BCs (calling subfunction
%     ApplyNaturalBC.m and ApplyEssBC.m). 
%     last update: 17 Nov 2015 H. Ritz; Y. Xu  
% Project Updates: Pulled from LinElast code base from MAE-5700 course,
%     updates were made to implement Kirchhoff Theory Plate Bending
% Update Authors: Dalton and Sairam
% *************************************************************************

function [globalSystem,boundStruct] = ApplyBC(boundStruct,meshStruct,globalSystem)
%% Apply natural BC(s)
% Loop over all natural boundaries
for i = 1:size(boundStruct.SurfNat,1)
    % update global K and F
    globalSystem = ApplyNaturalBC(i,boundStruct,meshStruct,globalSystem);
end


%% Apply essential BC(s)

essDOF = []; % TessDOF stores all the essential DOF
essVals = [];  % TessVals stores all the corresponding applied essential values

% Loop over all essential boundaries 
for surface_no = 1:size(boundStruct.SurfEssV,1)
    [cessDOF, cessVals]=ApplyEssBC(surface_no,boundStruct);% assign essential BCs
    essDOF =  [essDOF; cessDOF];
    essVals = [essVals; cessVals];

end

%% Captruing unique essDOF

% Because the corner nodes are at two boundaries, we need to find the
% unique numbers of debc (you dont want to apply the same EBC twice!)
[essDOF,m,~] = unique(essDOF);
essVals = essVals(m);

%% Package variables in structs
boundStruct.essDOF = essDOF;
boundStruct.ebcVals = essVals;
