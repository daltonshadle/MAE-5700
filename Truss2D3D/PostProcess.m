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
        error('Missing code. Please fill in postprocessing for 3D elements.');
end
strain=sum(operator.*d(gatherMat),2)./L;   % element strain
stress=elYM.*strain;  % element stress
force=elArea.*stress; % internal element force

% Package variables into the output structs
globalSystem.strain=strain;
globalSystem.stress=stress;
globalSystem.force =force;



PresentResults(globalSystem,meshStruct,boundStruct); % print out the results of the problem, including plots
PlotTensionCompression(globalSystem,meshStruct);
