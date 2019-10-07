function [ke,fe] = BeamElem(elmID,meshStruct)
% [localstiffnessmatrix, localforcevector] = BeamElem(elementnumber,meshStruct)
% generate the local stiffness matrix and local force vector from
% distributed loads for use with BEAM code.
% last edit: 5 August 2015 H. Ritz

% unpack necessary input
elCon  =meshStruct.elCon;
nCoords=meshStruct.nCoords;
elEI   =meshStruct.elEI;
elDistLoad=meshStruct.elDistLoad;



gn1 = elCon(elmID,1); % Extract the global node numbers
gn2 = elCon(elmID,2); % for the current element

L = nCoords(gn2)-nCoords(gn1);    % Extract the length of the element

EI = elEI(elmID);  % extract element property

% local stiffness matrix
ke = EI/L^3*[12,      6*L,     -12,      6*L;     
               6*L,    4*L^2,    -6*L,    2*L^2;
               -12,     -6*L,      12,     -6*L;
               6*L,    2*L^2,    -6*L,    4*L^2];

% assume constant distributed load and make the local force vector
p = elDistLoad(elmID);

fe = p*L/2*[1; L/6; 1; -L/6]; % Here, we assume constant q in each element
                              % Adjust if otherwise.  