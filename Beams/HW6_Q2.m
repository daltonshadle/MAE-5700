% *************************************************************************
% MAE-5700 HW6 Q2
% Authors: Sairam and Dalton
% Date: 10/17/19
% Notes: Lines of code listed here are the changes made to the original
%        Beam package, only the changes are listed.
% *************************************************************************

% From InputData.m ********************************************************
% ...
% follows the form [startLoad, endLoad, startPOI, endPOI]
spanVaryDistLoad = [0, 5000, 1, 3];
elVDL = zeros(size(nCoords));

% iterate over each varying distributed load (VDL)
for i=1:size(spanVaryDistLoad,1)
    % get the length of the VDL
    startPOI = pointsOfInterest(spanVaryDistLoad(i,3));
    endPOI = pointsOfInterest(spanVaryDistLoad(i,4));
    vdlLength = endPOI - startPOI;
    
    % get the starting load of the vdl (beta in our equation)
    vdlStart = spanVaryDistLoad(i,1);
    
    % get the vdl per unit length (alpha in our equation)
    vdlPerLength = (spanVaryDistLoad(i,2) - spanVaryDistLoad(i,1)) / vdlLength;
    
    % get all coordinates in this vdl
    vldCoordCond = (nCoords >= startPOI & nCoords <= endPOI);
    vdlCoords = nCoords(vldCoordCond);
    
    % find the vdl conditions
    tempElemVDL = vdlPerLength * vdlCoords + vdlStart;
    
    % assign to total;
    elVDL(vldCoordCond) = elVDL(vldCoordCond) + tempElemVDL;
end
% ...
meshStruct.elVDL = elVDL;


% From BeamElem.m *********************************************************
% ...
% unpack necessary input
elVDL =meshStruct.elVDL;
% ...

% assume constant distributed load and make the local force vector
p = elDistLoad(elmID);

fe_const = p*L/2*[1; L/6; 1; -L/6]; % Here, we assume constant q in each element
                                    % Adjust if otherwise.  
                                    
% assume varying distributed load and make the local force vector
q1 = elVDL(gn1);
q2 = elVDL(gn2);
fe_vdl = L/20*[7*q1 + 3*q2; 
               L*(q1 + (2/3)*q2); 
               3*q1 + 7*q2; 
               -L*((2/3)*q1 + q2)];

% create total fe
fe = fe_const + fe_vdl;