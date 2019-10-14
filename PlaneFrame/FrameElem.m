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

error('Generate the local stiffness matrix and local force vector, both in global coordinates')