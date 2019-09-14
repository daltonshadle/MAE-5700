function PlotTensionCompression(globalSystem,meshStruct)
% PLOTTENSIONCOMPRESSION(globalSystem,meshStruct)
% This function makes a plot of the truss using blue for elements in
% tension, red for elements in compression, and green for no-load members
% last edit: 21 November 2017 H. Ritz

% unpack necessary input
numEls    =meshStruct.numEls;
numDim    =meshStruct.numDim;
elCon     =meshStruct.elCon;
nCoords   =meshStruct.nCoords;
stress    =globalSystem.stress;

fighand=figure; % start a figure window and set some default properties
set(fighand,'defaultLineLineWidth',3)
set(fighand,'defaultTextFontSize',24)
set(fighand,'defaultAxesFontSize',24)
set(fighand,'defaultAxesFontWeight','bold')

% identify which elements have which sense of stress
eps=max(abs(stress))/10^14; % expected accuracy of the code
tensionEls=find(stress>eps);
compressionEls=find(stress<-eps);
noloadEls=find(abs(stress)<eps);

% make sure all elements are accounted for
if (length(tensionEls)+length(compressionEls)+length(noloadEls)~=numEls)
    display('Missing some elements?')
end

hold on;
switch numDim
    case 2
        % force a fixed legend
        plot(nan,nan,'-b',nan,nan,'-r',nan,nan,'-g');
        legend('Tension','Compression','No Load','AutoUpdate','off');
        for i = 1:numEls % loop over the elements
            gn1 = elCon(i,1);          % Extract the global node numbers
            gn2 = elCon(i,2);          % for the current element
            
            % get the initial coordinates and plot the truss
            XX = nCoords([gn1, gn2],1);
            YY = nCoords([gn1, gn2],2);
            if ismember(i,tensionEls)
                plotstyle='-b';
            elseif ismember(i,compressionEls)
                plotstyle='-r';
            elseif ismember(i,noloadEls)
                plotstyle='-g';
            else
                error('Something went wrong.')
            end
            plot(XX,YY,plotstyle);hold on;
        end
    case 3
        % force a fixed legend
        plot3(nan,nan,'-b',nan,nan,'-r',nan,nan,'-g');
        legend('Tension','Compression','No Load','AutoUpdate','off');
        for i = 1:numEls % loop over the elements
            gn1 = elCon(i,1);          % Extract the global node numbers
            gn2 = elCon(i,2);          % for the current element
            
            % get the initial coordinates and plot the truss
            XX = nCoords([gn1, gn2],1);
            YY = nCoords([gn1, gn2],2);
            ZZ = nCoords([gn1, gn2],3);
            if ismember(i,tensionEls)
                plotstyle='-b';
            elseif ismember(i,compressionEls)
                plotstyle='-r';
            elseif ismember(i,noloadEls)
                plotstyle='-g';
            else
                error('Something went wrong.')
            end
            plot3(XX,YY,ZZ,plotstyle);hold on;
        end
end
% label the plot
title('Stress Sense');
axis equal
