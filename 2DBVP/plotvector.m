function plotvector(arrows,points,title_string,scale);
% It plots the flux vector at Gauss points
% last update: 21 November 2017 H. Ritz


   figure; 
   fprintf(1, ' Plotting vector...\n');
   
   if nargin == 3
       scale = 0.9;
   end
    
   quiver(points(:,1),points(:,2), arrows(:,1), arrows(:,2),'LineWidth',1.0,'AutoScaleFactor',scale);
   hold off;
title(title_string); xlabel('x'); ylabel('y'); axis image;

