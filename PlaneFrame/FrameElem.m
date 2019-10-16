function [ke,fe] = FrameElem(elmID,meshStruct)
% [localstiffnessmatrix, localforcevector] = FrameElem(elementnumber,meshStruct)
% generate the local stiffness matrix and local force vector from
% distributed loads for use with FRAME code.
% last edit: 5 August 2015 H. Ritz

% unpack necessary input
elCon  =meshStruct.elCon;
nCoords=meshStruct.nCoords;
elEI   =meshStruct.elEI;
elEA   =meshStruct.elEA;
elDistLoad=meshStruct.elDistLoad;

gn1 = elCon(elmID,1); % Extract the global node numbers 
gn2 = elCon(elmID,2); % for the current element

x1 = nCoords(gn1,1); y1=nCoords(gn1,2); % extract the nodal coordinates for
x2 = nCoords(gn2,1); y2=nCoords(gn2,2); % each node in the current element

L = sqrt((x2-x1)^2 + (y2-y1)^2); % L = length of the element 
c = (x2-x1)/L; % cosine of the angle of the element with the X axis
s = (y2-y1)/L; % sine of the angle of the element with the X axis


% create rotation matrix to transform from local to global and visa-versa
R = [ c, s, 0, 0, 0, 0;
     -s, c, 0, 0, 0, 0;
      0, 0, 1, 0, 0, 0;
      0, 0, 0, c, s, 0;
      0, 0, 0,-s, c, 0;
      0, 0, 0, 0, 0, 1;];

% initialize truss and beam variables
EI = elEI(elmID);
Kt = elEA(elmID)/L;

% create truss elemental stiffness matrix in local coord
Ke_truss = [ Kt, 0, 0, -Kt, 0, 0;
              0, 0, 0,   0, 0, 0;
              0, 0, 0,   0, 0, 0;
            -Kt, 0, 0,  Kt, 0, 0;
              0, 0, 0,   0, 0, 0;
              0, 0, 0,   0, 0, 0;];
Ke_beam = EI/L^3*[ 0,   0,     0, 0,    0,     0;
                   0,  12,   6*L, 0,  -12,   6*L;     
                   0, 6*L, 4*L^2, 0, -6*L, 2*L^2;
                   0,   0,     0, 0,    0,     0;
                   0, -12,  -6*L, 0,   12,  -6*L;
                   0, 6*L, 2*L^2, 0, -6*L, 4*L^2];
Ke_local = Ke_truss + Ke_beam;

% transform to global coordinate system
ke = R' * Ke_local * R;

% assume constant distributed load and make the local force vector
p = elDistLoad(elmID);

% Assemble local elemental force, we assume constant q in each element
Fe_local = p*L/2*[0; 1; L/6; 0; 1; -L/6]; 

% Transform to global coord
fe = R' * Fe_local;
