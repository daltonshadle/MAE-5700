% From InputData.m
%...
elThermCoef = 0*ones(numEls,1); % Thermal Coefficient of Expansion
elDeltaTemp = 0*ones(numEls,1);    % Delta Temperature
%...
meshStruct.elThermCoef =elThermCoef;
meshStruct.elDeltaTemp =elDeltaTemp;


% From TrussElem.m
function [ke, fe_therm] = TrussElem(elmID,meshStruct)
%...
elThermCoef = meshStruct.elThermCoef;
elDeltaTemp = meshStruct.elDeltaTemp;
%...
case 2 % 2D problems        
    %...
    const_therm = elArea(elmID)*elYM(elmID)*elThermCoef(elmID)*elDeltaTemp(elmID); % constant coefficient for the element
    fe_therm = const_therm * [-c; -s; c; s]; % contribution of force vector from thermal effects
case 3 % 3D problems
    %...
    const_therm = elArea(elmID)*elYM(elmID)*elThermCoef(elmID)*elDeltaTemp(elmID); % constant coefficient for the element
    %...
    fe_therm = const_therm * [-a; -b; -c; a; b; c]; % contribution of force vector from thermal effects

    
% From Assembly.m
%...
f_total_therm = zeros(numNodes * numDOF, 1);
%...
[ke, fe_therm] = TrussElem(e,meshStruct); % make the local stiffness matrix
%...
local2global_map = gatherMat(e, :);
f_total_therm(local2global_map) = f_total_therm(local2global_map) + fe_therm;
%...
% Package variables into the output structs
globalSystem.f_total_therm = f_total_therm;


% From Soln.m
%...
f_total_therm = globalSystem.f_total_therm;
%...
f_Th    = f_total_therm(indF); % Extract force from thermal effects
%...
% solve for d_F
d_F	=K_F\( f_F + f_Th - K_EF'* d_E);
%...
% compute the reaction forces on the DOF with essential BCs
reactionVec = K_E*d_E+K_EF*d_F-F(essDOF) - f_total_therm(essDOF);


% From PostProcess.m
%...
elThermCoef = meshStruct.elThermCoef;
elDeltaTemp = meshStruct.elDeltaTemp;
%...
% Calculating all strains
strain=sum(operator.*d(gatherMat),2)./L;   % element strain
strain_therm = elThermCoef.*elDeltaTemp.*ones(numEls,1);
strain_mech = strain - strain_therm;

% Calculating all stresses
stress=elYM.*strain;  % element stress
stress_therm = elYM.*strain_therm;
stress_mech = stress - stress_therm;

% Calculating all forces
force=elArea.*stress; % internal element force
force_therm = elArea.*stress_therm;
force_mech = force - force_therm;

% Package variables into the output structs
globalSystem.strain =strain;
globalSystem.strain_mech =strain_mech;
globalSystem.strain_therm =strain_therm;
globalSystem.stress =stress;
globalSystem.stress_mech= stress_mech;
globalSystem.stress_therm =stress_therm;
globalSystem.force =force;
globalSystem.force_mech =force_mech;
globalSystem.force_therm =force_therm;

