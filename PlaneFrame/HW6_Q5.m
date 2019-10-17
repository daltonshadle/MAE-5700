% *************************************************************************
% MAE-5700 HW6 Q5
% Authors: Sairam and Dalton
% Date: 10/17/19
% Notes: Lines of code listed here are the changes made to the original
%        PlaneFrame package, only the changes are listed.
% *************************************************************************

% From FrameMesh.m ********************************************************
% ...
pointsOfInterest=[0 0 ;
                  0 3; 
                  0 6;
                  4 6;
                  4 0;
                  4 3];
% ...
spanCon=[1 2;
         2 3;
         2 6;
         3 4;
         4 6;
         6 5];
% ...
spanEls=[20 20 20 20 20 20];
% ...
maxLength=4;


% From InputData.m ********************************************************
% ...
spanDistributedLoad = [0 0 -2e3 -2e3 0 0]';
% ...
essBCs=[1 1 0;
        1 2 0;
        1 3 0;
        5 1 0;
        5 2 0;
        5 3 0]; 
% ...
appForces=[2 1 4000;
           3 1 4000];
% ...
for frc=1:size(appForces,1)
    % ...
    % find global node number, Note: Sairam wanted to note that this only 
    % took one line of code, not two lines
    gnn = find(poi(1,1)== nCoords(:,1) & poi(1,2) == nCoords(:,2));
    % ...
end


% From FrameElem.m ********************************************************
% ...
gn1 = elCon(elmID,1); % Extract the global node numbers
gn2 = elCon(elmID,2); % for the current element
x1=nCoords(gn1,1); y1=nCoords(gn1,2); % extract the nodal coordinates for
x2=nCoords(gn2,1); y2=nCoords(gn2,2); % each node in the current element
L=sqrt((x2-x1)^2+(y2-y1)^2); % L = length of the element
c=(x2-x1)/L; % cosine of the angle of the element with the X axis
s=(y2-y1)/L; % cosine of the angle of the element with the Y axis

EI = elEI(elmID); % extract element property EI
EA = elEA(elmID); % extract element property EA

% local stiffness matrix  in local co-ordinates
% ke1 = truss stiffness matrix, ke2 = beam stiffness matrix
ke1 = EA/L*[1 0 0 -1 0 0;0 0 0 0 0 0;0 0 0 0 0 0;-1 0 0 1 0 0; 0 0 0 0 0 0;0 0 0 0 0 0];
ke2 = EI/L^3*[0 0 0 0 0 0;0 12 6*L 0 -12 6*L;0 6*L 4*L^2 0 -6*L 2*L^2;     
              0 0 0 0 0 0;0 -12 -6*L 0 12 -6*L;0 6*L 2*L^2 0 -6*L 4*L^2];

% add beam and truss stiffness matrices to get total stiffness matrix
ke_local = ke1 + ke2;    

% initialize rotation matrix
Re = [c s 0 0 0 0;-s c 0 0 0 0;0 0 1 0 0 0;0 0 0 c s 0;0 0 0 -s c 0;0 0 0 0 0 1];

% transform stiffness matrix from local to global coord
ke = Re'*ke_local*Re;

% assume constant distributed load and make the local force vector
p = elDistLoad(elmID);
fe_local = p*L/2*[0; 1; L/6; 0; 1; -L/6];

% transform elemental force vector from local to global coord
fe = Re'*fe_local;



% From Soln.m *************************************************************
% ...
for ebc=1:numEBC
    % ...
    % find global node number, Note: Sairam wanted to note that this only 
    % took one line of code, not two lines
    gnn = find(poi(1,1)== nCoords(:,1) & poi(1,2) == nCoords(:,2));
    % ...
end


% From PostProcess.m ******************************************************
% ...
% initialize variables for plotting and results
deflection_x = [];
deflection_y = [];
axialForce = zeros(numEls,1);
shearForce = zeros(numEls,1);
bendingMoment1 = zeros(numEls,1);
bendingMoment2 = zeros(numEls,1);
% ...
for elmID = 1:numEls      % loop over elements to postprocess and plot
    % ...
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
    allXi = linspace(-1,1,nplot);
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
    end
    deflection_x    = [deflection_x;deflections(1,:)]; 
    deflection_y   = [deflection_y;deflections(2,:)]; 
    axialForce(elmID) = elEA(elmID)/L*[-1 0 0 1 0 0]*Re*elmSoln;
    shearForce(elmID) = -elEI(elmID)/J^3*[0 3/2 3*L/4 0 -3/2 3*L/4]*Re*elmSoln;
    bendingMoment1(elmID) = elEI(elmID)/J^2*[0 -3/2 -L 0 3/2 -L/2]*Re*elmSoln;
    bendingMoment2(elmID) = elEI(elmID)/J^2*[0 3/2 L/2 0 -3/2 L]*Re*elmSoln;
end
% ...
% package necessary variables
globalSystem.deflection_x   =deflection_x;
globalSystem.deflection_y   =deflection_y;
globalSystem.bendingMoment1 =bendingMoment1;
globalSystem.bendingMoment2 =bendingMoment2;
% ...
% print results
for i=1:numEls
    fprintf(FID,'%d\t%e\t%e\t%e\t%e\n',i,axialForce(i),shearForce(i),bendingMoment1(i),bendingMoment2(i));
end




