function [ke, fe_therm] = TrussElem(elmID,meshStruct)
% localstiffnessmatrix = TrussElem(elementnumber,meshStruct)
% generate the local stiffness matrix for use with TRUSS2D3D code.
% last edit: 30 July 2015 H. Ritz

% unpack necessary input
elCon  =meshStruct.elCon;
nCoords=meshStruct.nCoords;
elArea =meshStruct.elArea;
elYM   =meshStruct.elYM;
elThermCoef = meshStruct.elThermCoef;
elDeltaTemp = meshStruct.elDeltaTemp;

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
        const_therm = elArea(elmID)*elYM(elmID)*elThermCoef(elmID)*elDeltaTemp(elmID); % constant coefficient for the element
        ke = ...
            const*[c2   c*s   -c2   -c*s;  % stiffness matrix for 2D
            c*s     s2     -c*s     -s2;   % truss element in the global
            -c2     -c*s     c2       c*s; % coordinate system
            -c*s    -s2      c*s      s2];
        fe_therm = const_therm * [-c; -s; c; s]; % contribution of force vector from thermal effects
    case 3 % 3D problems
        z1=nCoords(gn1,3); z2=nCoords(gn2,3); % extract z-coord for both nodes
        L=sqrt((x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2); % L = length of the element
        
        a = (x2-x1) / L; % direction cosine of projection on x-axis
        b = (y2-y1) / L; % direction cosine of projection on y-axis
        c = (z2-z1) / L; % direction cosine of projection on z-axis
        
        const = elArea(elmID)*elYM(elmID)/L; % constant stiffness coefficient for the element
        const_therm = elArea(elmID)*elYM(elmID)*elThermCoef(elmID)*elDeltaTemp(elmID); % constant coefficient for the element
        re = [a, b, c, 0, 0, 0;  % rotation matrix for 3D case
              0, 0, 0, a, b, c];
        
        ke = const * re' * [1, -1; -1, 1] * re; % elemental stiffness matrix 
                                                % for 3D truss system in
                                                % global coord system
        fe_therm = const_therm * [-a; -b; -c; a; b; c]; % contribution of force vector from thermal effects
end
