function ke = TrussElem(elmID,meshStruct)
% localstiffnessmatrix = TrussElem(elementnumber,meshStruct)
% generate the local stiffness matrix for use with TRUSS2D3D code.
% last edit: 30 July 2015 H. Ritz

% unpack necessary input
elCon  =meshStruct.elCon;
nCoords=meshStruct.nCoords;
elArea =meshStruct.elArea;
elYM   =meshStruct.elYM;
numDim =meshStruct.numDim;
nnpe   =meshStruct.nnpe;


gn1 = elCon(elmID,1); % Extract the global node numbers 
gn2 = elCon(elmID,2); % for the current element

x1=nCoords(gn1,1); y1=nCoords(gn1,2); % extract the nodal coordinates for
x2=nCoords(gn2,1); y2=nCoords(gn2,2); % each node in the current element
          
switch numDim % the local stiffness matrix has a different form
              % depending on the number of spatial dimensions
    case 2 % 2D problems        
        L=sqrt((x2-x1)^2+(y2-y1)^2); % L = length of the element
        c=(x2-x1)/L; % cosine of the angle of the element with the X axis
        c2 = c^2;
        s=(y2-y1)/L; % cosine of the angle of the element with the Y axis
        s2 = s^2;
        
        const = elArea(elmID)*elYM(elmID)/L; % constant coefficient for the element
        ke = ...
            const*[c2   c*s   -c2   -c*s;  % stiffness matrix for 2D
            c*s     s2     -c*s     -s2;   % truss element in the global
            -c2     -c*s     c2       c*s; % coordinate system
            -c*s    -s2      c*s      s2];
    case 3 % 3D problems
        error('Missing code. Please fill in elemental stiffness matrix for 3D elements.');
end
