function [u du] = Exact ( x )

% This function specifies the analytical solution u(x)
% and its derivative du/dx at point x (if they exist
% and are known). This is problem specific.

% Remember to use vectorization, e.g. use .*, ./ and .^
% when computing this solution.

u    =      sin ( 4 * pi * x );
du   =     4*pi*cos( 4 * pi * x );


% % TAPERED BAR EXAMPLE
% for ii=1:length(x)
%     if x(ii)<5 % left
%         du(ii)=(36-4.*x(ii))./(8*x(ii));
%         u(ii)=4.5*log(x(ii))-(1/2)*x(ii)-4.5*log(2)+(1/2)*2;
%     elseif x(ii)>=5 % right
%         du(ii)=(24-4*x(ii))/(8*x(ii));
%         u(ii)=3*log(x(ii))-(1/2)*x(ii)+1.5*log(5)-4.5*log(2)+1;
%     end
% end

