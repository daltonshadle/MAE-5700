% ApplyBC
% This function applies natural and essential BCs (calling subfunction
% ApplyNaturalBC.m and ApplyEssBC.m). 

% last update: 17 Nov 2015 H. Ritz; Y. Xu  

function [globalSystem,boundStruct] = ApplyBC(boundStruct,meshStruct,globalSystem)

%% Apply natural BC(s)
for i = 1:size(boundStruct.SurfNat,1)    % Loop over all natural boundaries
    globalSystem = ApplyNaturalBC(i,boundStruct,meshStruct,globalSystem);  % update global K and F

end

%% Apply essential BC(s)

essDOF = []; % TessDOF stores all the essential DOF
essVals = [];  % TessVals stores all the corresponding applied essential values

for surface_no = 1:size(boundStruct.SurfEssV,1)   % Loop over all essential boundaries 

    [cessDOF, cessVals]=ApplyEssBC(surface_no,boundStruct);% assign essential BCs
    essDOF =  [essDOF; cessDOF];
    essVals = [essVals; cessVals];

end
%% Captruing unique essDOF

% Because the corner nodes are at two boundaries, we need to find the
% unique numbers of debc (you dont want to apply the same EBC twice!)
[essDOF,m,~] = unique(essDOF);
essVals = essVals(m);

boundStruct.essDOF = essDOF;
boundStruct.ebcVals = essVals;
%%

