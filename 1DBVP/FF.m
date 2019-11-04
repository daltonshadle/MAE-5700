function value = FF ( x )
% Value = FF ( X )
% It returns the Value of the right hand side 
% in the ODE: -d/dx (p(x) du/dx)  =  f(x)
% You must change this for every problem you want to solve.
% last edit: 22 March 2015 H. Ritz

value = 16*pi^2*sin(4*pi*x);

% % Tapered bar problem
% value=8;