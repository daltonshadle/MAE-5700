function globalSystem=PostProcess(globalSystem,meshStruct,boundStruct);
% globalSystem=POSTPROCESS(globalSystem,meshStruct,boundStruct)
% Calculate bending moment, shear force, and axial force for FRAME problem.
% Plot the deformed shapes of the frame.
% Print relevant problem data and results.
% last edit: 5 August 2015 H. Ritz


% unpack necessary input
numEls    =meshStruct.numEls;
gatherMat =meshStruct.gatherMat;
elEI      =meshStruct.elEI;
elEA      =meshStruct.elEA;
nCoords   =meshStruct.nCoords;
elCon     =meshStruct.elCon;
d         =globalSystem.d;

numSpans =meshStruct.numSpans;
numEq    =meshStruct.numEq;
numNodes =meshStruct.numNodes;
numDOF    =meshStruct.numDOF;

appForces=boundStruct.appForces;
essBCs   =boundStruct.essBCs;
reactionVec=globalSystem.reactionVec;

% want to plot both the initial and deformed shapes
fighand=figure; % start a figure window and set some default properties
set(fighand,'defaultLineLineWidth',3)
set(fighand,'defaultTextFontSize',24)
set(fighand,'defaultAxesFontSize',24)
set(fighand,'defaultAxesFontWeight','bold')
% need a magnification factor for deformed shape since
% deformations are so small. You may need to adjust this for different
% types of problems.
tmpd=[reshape(d,3,[])]';
maxd=max(max(abs(tmpd(:,1:2))));
if maxd==0
    magFac=1;
else
    magFac=0.1*max(max(nCoords))/maxd;
end
deflection_x = [];
deflection_y = [];
axialForce = zeros(numEls,1);
shearForce = zeros(numEls,1);
bendingMoment1 = zeros(numEls,1);
bendingMoment2 = zeros(numEls,1);
for elmID = 1:numEls      % loop over elements to postprocess and plot
    % To find the deformation of each point in the frame we need to find
    % the axial deformation and the transverse deflection separately from
    % each other. The methodology should be to find several points, xi,
    % along the length of each element and use the shape functions from the
    % beam FE model calculated at those xi to find the transverse
    % deflections in the local coordinates. Axial deflection should vary
    % linearly from one endpoint to the other, so use the shape functions
    % from the truss elements discussed previously. Use rotMat' to
    % transform those deflections into global coordinates.
    elmSoln=d(gatherMat(elmID,:)); % extract element nodal displacements
    EI   = elEI(elmID);            % Young's modulus x I
    gn1 = elCon(elmID,1); % Extract the global node numbers
    gn2 = elCon(elmID,2); % for the current element
    x1=nCoords(gn1,1);
    y1=nCoords(gn1,2); % extract the nodal coordinates for
    x2=nCoords(gn2,1);
    y2=nCoords(gn2,2); % each node in the current element
    L=sqrt((x2-x1)^2+(y2-y1)^2); % L = length of the element
    c=(x2-x1)/L; % cosine of the angle of the element with the X axis
    s=(y2-y1)/L; % cosine of the angle of the element with the Y axis
    Re = [c s 0 0 0 0;-s c 0 0 0 0;0 0 1 0 0 0;0 0 0 c s 0;0 0 0 -s c 0;0 0 0 0 0 1];
    J  = L / 2;         % compute Jacobian for transformation.
    nplot=3; %number of points along each element for plotting
    % Compute displacements, moments and shear forces
    allXi = linspace(-1,1,nplot)
    for ind=1:nplot
        xi=allXi(ind); % for each xi ...
        xy(1,ind) = (1-xi)/2*x1 + (1+xi)/2*x2;
        xy(2,ind) = (1-xi)/2*y1 + (1+xi)/2*y2;       
        % shape functions evaluated at this point
        N1 = [(1-xi)/2, 0, 0, (1+xi)/2, 0 , 0];
        N2 = [0, 1/4*(1-xi)^2*(2+xi), L/8*(1-xi)^2*(1+xi), 0, ...
            1/4*(1+xi)^2*(2-xi), L/8*(1+xi)^2*(xi-1)];
        deflections(1,ind) = N1*elmSoln;
        deflections(2,ind) = N2*elmSoln;
        
        %             'xy should have dimensions 2 X nplot.\n',...
        %             '\txy(1,ind)=global x cocordinate of this xi location.\n',...
        %             '\txy(2,ind)=global y cocordinate of this xi location.\n',...
        %             'deflections should have dimensions 2 X nplot.\n',...
        %             '\tdeflections(1,ind)=ux in global coordinates of this xi location.\n',...
        %             '\tdeflections(2,ind)=uy in global coordinates of this xi location.\n']))
    end
    deflection_x    = [deflection_x;deflections(1,:)]; 
    deflection_y   = [deflection_y;deflections(2,:)]; 
    axialForce(elmID) = elEA(elmID)/L*[-1 0 0 1 0 0]*Re*elmSoln;
    shearForce(elmID) = elEI(elmID)/J^3*[0 3/2 3*L/4 0 -3/2 3*L/4]*Re*elmSoln;
    bendingMoment1(elmID) = elEI(elmID)/J^2*[0 -3/2 -L 0 3/2 L/2]*Re*elmSoln;
    bendingMoment2(elmID) = elEI(elmID)/J^2*[0 3/2 L/2 0 -3/2 L]*Re*elmSoln;
    % plot the undeformed frame
    plot([x1 x2],[y1 y2],'b-');hold on;
    xy(1,:)
    xy(2,:)
    % plot the deformed frame, using the magnifaction factor when
    % adding the displacements to the original positions.
    plot(xy(1,:)+magFac*deflections(1,:),xy(2,:)+magFac*deflections(2,:),'r-');hold on;
end
axis equal
% label the plot
title('Frame Plot');
legend('Initial',['Deformed (',num2str(magFac),'X)'])
% Package variables into the output struct
globalSystem.deflection_x   =deflection_x;
globalSystem.deflection_y   =deflection_y;
globalSystem.bendingMoment1 =bendingMoment1;
globalSystem.bendingMoment2 =bendingMoment2;
globalSystem.shearForce    =shearForce;
globalSystem.axialForce    =axialForce;

% Print problem input and results
FID=1; % FID=1 prints to the screen.
% FID=fopen('FEAoutput.txt','w'); % Use this line if you want to
% print to a file using FOPEN


% print mesh parameters
fprintf(FID,'\n\tFrame Parameters \n\t----------------\n');
fprintf(FID,'No. of Spans     %d \n',numSpans);
fprintf(FID,'No. of Elements  %d \n',numEls);
fprintf(FID,'No. of Nodes     %d \n',numNodes);
fprintf(FID,'No. of Equations %d \n',numEq);

% find the directions of the applied forces
if ~isempty(appForces)
    appdir(appForces(:,2)==1)='x';
    appdir(appForces(:,2)==2)='y';
    appdir(appForces(:,2)==3)='M';
    % print the applied forces
    fprintf(FID,'\n\n\tApplied Forces \n\t--------------\n');
    fprintf(FID,'node #\tdir. of app. force\tvalue of app. force\n');
    for i=1:size(appForces,1)
        fprintf(FID,'%d\t\t%s\t\t%e\n',appForces(i,1),appdir(i),appForces(i,3));
    end
end
% print the solution vector
fprintf(FID,'\n\n\tNodal Displacements \n\t-------------------\n');
fprintf(FID,'node #\tx-displacement\ty-displacement\trotation\n');
fprintf(FID,'%d\t%e\t%e\t%e\n',[1:numNodes;reshape(d,numDOF,numNodes)]);

% find the directions of the essential BC
rxndir(essBCs(:,2)==1)='x';
rxndir(essBCs(:,2)==2)='y';
rxndir(essBCs(:,2)==3)='M';
% print the reaction forces
fprintf(FID,'\n\n\tReaction Forces\n\t---------------\n');
fprintf(FID,'POI #\tdir. of rxn\tvalue of rxn\n');
for i=1:size(essBCs,1)
    fprintf(FID,'%d\t%s\t\t%e\n',essBCs(i,1),rxndir(i),reactionVec(i));
end

% print the strains, stresses, and internal forces for each element
fprintf(FID,'\n\n\tElement Results \n\t---------------\n');
fprintf(FID,'el #\taxial\t\tshear\t\tmoment 1\tmoment 2\n');
%fprintf(FID,'%d\t%e\t%e\t%e\t%e\n',[1:numEls;axialForce;shearForce;bendingMoment1;bendingMoment2]);
for i=1:numEls
    fprintf(FID,'%d\t%e\t%e\t%e\t%e\n',i,axialForce(i),shearForce(i),bendingMoment1(i),bendingMoment2(i));
end

