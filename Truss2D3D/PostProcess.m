function globalSystem = PostProcess(globalSystem,meshStruct,boundStruct)
% globalSystem = POSTPROCESS(globalSystem,meshStruct,boundStruct)
% Calculate strains, stresses, and internal forces for each element in the
% TRUSS2D3D problem. Plot the intial and deformed shapes of the truss. 
% Print relevant problem data and results.
% last edit: 30 July 2015 H. Ritz

% unpack necessary input
numEls    =meshStruct.numEls;
numDim    =meshStruct.numDim;
elCon     =meshStruct.elCon;
gatherMat =meshStruct.gatherMat;
nCoords   =meshStruct.nCoords;
elYM      =meshStruct.elYM;
elArea    =meshStruct.elArea;
elThermCoef =meshStruct.elThermCoef ;
elDeltaTemp = meshStruct.elDeltaTemp;

d         =globalSystem.d';

% initialize output vectors
strain=zeros(numEls,1);
stress=zeros(numEls,1);
force =zeros(numEls,1);

gn1 = elCon(:,1);  % Extract the global node numbers
gn2 = elCon(:,2);  % for all elements

switch numDim % post process differently, depending on spatial dimensions
    case 2 % 2D trusses
        % get nodal coordinates
        x1=nCoords(gn1,1); y1=nCoords(gn1,2);
        x2=nCoords(gn2,1); y2=nCoords(gn2,2);
        
        L=sqrt((x2-x1).^2+(y2-y1).^2); % L = initial length of the elements
        c=(x2-x1)./L; % cosine of the angle of the elements with the X axis
        s=(y2-y1)./L; % cosine of the angle of the elements with the Y axis
        operator=[-c -s c s];
    case 3 %3D trusses
        x1=nCoords(gn1,1); x2=nCoords(gn2,1); % extract x-coord for both nodes
        y1=nCoords(gn1,2); y2=nCoords(gn2,2); % extract y-coord for both nodes
        z1=nCoords(gn1,3); z2=nCoords(gn2,3); % extract z-coord for both nodes
        
        L=sqrt((x2-x1).^2 + (y2-y1).^2 + (z2-z1).^2); % L = length of the element
        
        a = (x2-x1) / L; % direction cosine of projection on x-axis
        b = (y2-y1) / L; % direction cosine of projection on y-axis
        c = (z2-z1) / L; % direction cosine of projection on z-axis
        
        operator = [-a, -b, -c, a, b, c];
end

% Calculating all strains
disp(size(operator));
disp(size(d(gatherMat)));
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


PresentResults(globalSystem,meshStruct,boundStruct); % print out the results of the problem, including plots
PlotTensionCompression(globalSystem,meshStruct);
