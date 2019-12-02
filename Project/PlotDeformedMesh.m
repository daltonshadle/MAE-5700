% ------------------------------------------------------------------------|
%                                                                         |
% MAE4700-5700, Finite Element Analysis for Mechanical & Aerospace Design |
%                                                                         |
% Copyright: Cornell University (this software should not be used without |
% written permission)                                                     |
%                                                                         |
% Authors: N. Zabaras (zabaras@cornell.edu) & Xiang Ma (xm25@cornell.edu) |
%                                                                         |
% ------------------------------------------------------------------------|
%
 function PlotDeformedMesh(string_title,MeshStruct,d)

nno=MeshStruct.numNodes;
nsd=MeshStruct.nsd;
Nodes=MeshStruct.nCoords';
Elems=MeshStruct.elCon';
nel=MeshStruct.numEls;
% After computing the displacements, we plot
% the given mesh in its initial and deformed shapes. 

displacements=[reshape(d,nsd,[])]';

    a = zeros(nsd);     % Find the scale factor for deformation shape
    for i=1:nsd
        a(i) = max(Nodes(i,:)) - min(Nodes(i,:));
    end
    aa = max(a);
    b = zeros(nsd);
    scale = 0;

    for i=1:nsd
        b(i) = max(displacements(:,i)) - min(displacements(:,i));
    end
    bb = max(b);
    scale = 0.1 * aa / bb;

    figure;
    hold on;

    for i= 1:nel

        glb = Elems(:,i)';
        XX = [Nodes(1, glb),Nodes(1,glb(1))];
        YY = [Nodes(2, glb),Nodes(2,glb(1))];
        plot(XX,YY,'b');hold on;

        displacedXX = Nodes(1, glb) + scale* displacements(glb,1)';
        displacedYY = Nodes(2, glb) + scale* displacements(glb,2)';
        displacedXX = [displacedXX, displacedXX(1)];
        displacedYY = [displacedYY, displacedYY(1)];
        plot (displacedXX,displacedYY,'LineStyle','--','Color','r');

    end

    legend('Initial shape','Deformed shape'); axis image;
    title(string_title); xlabel('x'); ylabel('y');

