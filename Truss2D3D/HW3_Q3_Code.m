% FROM TrussElem.m
case 3 % 3D problems
        z1=nCoords(gn1,3); z2=nCoords(gn2,3); % extract z-coord for both nodes
        L=sqrt((x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2); % L = length of the element
        
        a = (x2-x1) / L; % direction cosine of projection on x-axis
        b = (y2-y1) / L; % direction cosine of projection on y-axis
        c = (z2-z1) / L; % direction cosine of projection on z-axis
        
        const = elArea(elmID)*elYM(elmID)/L; % constant stiffness coefficient for the elemen
        
        re = [a, b, c, 0, 0, 0;  % rotation matrix for 3D case
              0, 0, 0, a, b, c];
        
        ke = const * re' * [1, -1; -1, 1] * re; % elemental stiffness matrix 
                                                % for 3D truss system in
                                                % global coord system

% From PostProcess.m
case 3 %3D trusses
        x1=nCoords(gn1,1); x2=nCoords(gn2,1); % extract x-coord for both nodes
        y1=nCoords(gn1,2); y2=nCoords(gn2,2); % extract y-coord for both nodes
        z1=nCoords(gn1,3); z2=nCoords(gn2,3); % extract z-coord for both nodes
        
        L=sqrt((x2-x1).^2 + (y2-y1).^2 + (z2-z1).^2); % L = length of the element
        
        a = (x2-x1) ./ L; % direction cosine of projection on x-axis
        b = (y2-y1) ./ L; % direction cosine of projection on y-axis
        c = (z2-z1) ./ L; % direction cosine of projection on z-axis
        
        operator = [-a, -b, -c, a, b, c];