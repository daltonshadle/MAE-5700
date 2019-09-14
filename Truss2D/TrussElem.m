function ke = TrussElem(elmID,meshStruct)
% localstiffnessmatrix = TRUSSELEM(elementnumber,meshStruct)
% generate the local stiffness matrix for use with TRUSS2D code.
% NOTE! as written this only generates the local stiffness matrix for 2D
% truss elements with 2 nodes per element and 2 DOF per node. A
% SWITCH/CASE structure could be used to incorporate local stiffness
% matrices for other types of truss elements
% last edit: 10 July 2015 H. Ritz

% unpack necessary input
elCon  =meshStruct.elCon;
nCoords=meshStruct.nCoords;
elArea =meshStruct.elArea;
elYM   =meshStruct.elYM;


gn1 = elCon(elmID,1); % Extract the global node numbers 
gn2 = elCon(elmID,2); % for the current element

x1=nCoords(gn1,1); y1=nCoords(gn1,2); % extract the nodal coordinates for
x2=nCoords(gn2,1); y2=nCoords(gn2,2); % each node in the current element

L=sqrt((x2-x1)^2+(y2-y1)^2); % L = length of the element 
c=(x2-x1)/L; % cosine of the angle of the element with the X axis
c2 = c^2;
s=(y2-y1)/L; % sine of the angle of the element with the X axis
s2 = s^2;           
    
const = elArea(elmID)*elYM(elmID)/L; % constant coefficient for the element
   ke = const*[c2      c*s    -c2     -c*s;  % stiffness matrix for 2D 
               c*s     s2     -c*s     -s2;  % truss element in the global
              -c2     -c*s     c2       c*s; % coordinate system
              -c*s    -s2      c*s      s2];
           