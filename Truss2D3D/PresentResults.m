function PresentResults(globalSystem,meshStruct,boundStruct)
% PresentResults(globalSystem,meshStruct,boundStruct)
% Print relevant problem data and results for TRUSS2D.
% last edit: 30 July 2015 H. Ritz

% unpack necessary input
numEls    =meshStruct.numEls;
numNodes  =meshStruct.numNodes;
numDim    =meshStruct.numDim;
numDOF    =meshStruct.numDOF;
numEq     =meshStruct.numEq;
elCon     =meshStruct.elCon;
gatherMat =meshStruct.gatherMat;
nCoords   =meshStruct.nCoords;
d         =globalSystem.d;
reactionVec= globalSystem.reactionVec; 
strain    =globalSystem.strain;      
stress    =globalSystem.stress;     
force     =globalSystem.force;  
appForces =boundStruct.appForces;
essBCs    =boundStruct.essBCs;

% Plot both the initial and deformed shapes
fighand=figure; % start a figure window and set some default properties
set(fighand,'defaultLineLineWidth',3)
set(fighand,'defaultTextFontSize',24)
set(fighand,'defaultAxesFontSize',24)
set(fighand,'defaultAxesFontWeight','bold')
% need a magnification factor for deformed shape since 
% deformations are so small. You may need to adjust this for different
% types of problems. 
if max(abs(d)) ~= 0
    magFac=0.1*max(max(abs(nCoords)))/max(abs(d));
else
    magFac=1;
end
gn1 = elCon(:,1); % Extract the global node numbers
gn2 = elCon(:,2);

switch numDim
    case 2
        for i = 1:numEls % loop over the elements
            % get the initial coordinates and plot the truss
            XX = nCoords([gn1(i), gn2(i)],1);
            YY = nCoords([gn1(i), gn2(i)],2);
            plot(XX,YY,'b-');hold on;
            % label the node numbers
            text(XX(1),YY(1),sprintf('%0.5g',gn1(i)));
            text(XX(2),YY(2),sprintf('%0.5g',gn2(i)));
            % get the deformed coordinates and plot the deformed shape
            disps=reshape(d(gatherMat(i,:)),2,2); % this has the displacements of the
            % ends of the current element
            dx=disps(1,:)'; % x displacements
            dy=disps(2,:)'; % y displacements
            % plot, using the magnifaction factor when adding the
            % displacements to the original nodal positions.
            plot(XX+magFac*dx,YY+magFac*dy,'r--');hold on;
        end
    case 3
        for i = 1:numEls % loop over the elements
            % get the initial coordinates and plot the truss
            XX = nCoords([gn1(i), gn2(i)],1);
            YY = nCoords([gn1(i), gn2(i)],2);
            ZZ = nCoords([gn1(i), gn2(i)],3);
            plot3(XX,YY,ZZ,'b-');hold on;
            % label the node numbers
            text(XX(1),YY(1),ZZ(1),sprintf('%0.5g',gn1(i)));
            text(XX(2),YY(2),ZZ(2),sprintf('%0.5g',gn2(i)));
            % get the deformed coordinates and plot the deformed shape
            disps=reshape(d(gatherMat(i,:)),3,2); % this has the displacements of the
            % ends of the current element
            dx=disps(1,:)'; % x displacements
            dy=disps(2,:)'; % y displacements
            dz=disps(3,:)'; % z displacements
            % plot, using the magnifaction factor when adding the
            % displacements to the original nodal positions.
            plot3(XX+magFac*dx,YY+magFac*dy,ZZ+magFac*dz,'r--');hold on;
        end
end
axis equal
% label the plot
title('Truss Plot');
legend('Initial',['Deformed (',num2str(magFac),'X)'])
% Print problem input and results
FID=1; % FID=1 prints to the screen. 
% FID=fopen('Output.txt','w'); % use this code if you want to print to
%                                a file using FOPEN, instead

% print mesh parameters
fprintf(FID,'\n\tTruss Parameters \n\t----------------\n');
fprintf(FID,'No. of Elements  %d \n',numEls);
fprintf(FID,'No. of Nodes     %d \n',numNodes);
fprintf(FID,'No. of Equations %d \n',numEq);

if ~isempty(appForces)
    % find the directions of the applied forces
    appdir(appForces(:,2)==1)='x';
    appdir(appForces(:,2)==2)='y';
    appdir(appForces(:,2)==3)='z';
    % print the applied forces
    fprintf(FID,'\n\n\tApplied Forces \n\t--------------\n');
    fprintf(FID,'node #\tdir. of app. force\tvalue of app. force\n');
    for i=1:size(appForces,1)
        fprintf(FID,'%d\t\t%s\t\t%e\n',appForces(i,1),appdir(i),appForces(i,3));
    end
end

% print the solution vector
fprintf(FID,'\n\n\tNodal Displacements \n\t-------------------\n');
nodaldisplacements=[1:numNodes;reshape(d,numDOF,numNodes)];
switch numDim
    case 2
        fprintf(FID,'node #\tx-displacement\ty-displacement\n');
        fprintf(FID,'%d\t%e\t%e\n',nodaldisplacements);
    case 3
        fprintf(FID,'node #\tx-displacement\ty-displacement\tz-displacement\n');
        fprintf(FID,'%d\t%e\t%e\t%e\n',nodaldisplacements);
end

% find the directions of the essential BC
rxndir(essBCs(:,2)==1)='x';
rxndir(essBCs(:,2)==2)='y';
rxndir(essBCs(:,2)==3)='z';
% print the reaction forces
fprintf(FID,'\n\n\tReaction Forces \n\t---------------\n');
fprintf(FID,'node #\tdir. of rxn\tvalue of rxn\n');
for i=1:size(essBCs,1)
    fprintf(FID,'%d\t%s\t\t%e\n',essBCs(i,1),rxndir(i),reactionVec(i));
end

% print the strains, stresses, and internal forces for each element
fprintf(FID,'\n\n\tElement Results \n\t---------------\n');
fprintf(FID,'el #\tstress\t\tstrain\t\tinternal force\n');
fprintf(FID,'%d\t%e\t%e\t%e\n',[1:numEls;stress';strain';force']);

