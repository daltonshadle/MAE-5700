% ------------------------------------------------------------------------|
%                                                                         |
% MAE4700-5700, Finite Element Analysis for Mechanical & Aerospace Design |
%                                                                         |
% Copyright: Cornell University (this software should not be used without |
% written permission)                                                     |
%                                                                         |
% Authors: N. Zabaras (zabaras@cornell.edu) & Xiang Ma (xm25@cornell.edu) |
% Last update: 9 November 2015 Y. Xu                                      |                                  |
% ------------------------------------------------------------------------|
%
function [globalSystem,boundStruct] = ApplyBC(boundStruct,meshStruct,globalSystem)
% This function applies both natural and essential BCs.

% unpack the things you need
if isempty(boundStruct.SurfEssV)
    SideEID = [];
else
    SideEID = boundStruct.SurfEssV(:,1); % surfaces that have essential BC
end
if isempty(boundStruct.SurfNat)
    SideNID = [];
else
    SideNID = boundStruct.SurfNat(:,1); % surfaces that have natural BC
end
% Apply natural BC
for i = 1:length(SideNID)    % Loop over all boundaries

    globalSystem = ApplyNaturalBC(i,boundStruct,meshStruct,globalSystem); 

end

% Apply essential BC
tdebc    = [];  tebcVals = [];     % Initialize tdebc, tebcVals

for i = 1:length(SideEID)    % Loop over all boundaries

    [bc,vals]=ApplyEssBC(i,boundStruct);% assign essential BCs
    tdebc =  [tdebc, bc];
    tebcVals = [tebcVals; vals];

end

% Because the corner nodes are at two boundaries, we need to find the
% unique numbers of debc (you dont want to apply the same EBC twice!)
[debc,m,~] = unique(tdebc);
ebcVals = tebcVals(m);
boundStruct.essDOF = debc;
boundStruct.ebcVals = ebcVals;
clear tdebc tebcVals;



