function [ex,S,M,w]=exact

% Analytical solution for the example problem 
% This function needs to change from problem to 
% problem - if you dont have analytical results, 
% remove all of these references to the exact solution. 
%
% last modified: 18 February 2015 H. Ritz

L=12; a1=4; p1=10; a2=8; p2=-5; a3=8; m=-20; p4=20; p3=1;
E=1e4; I=1;
i=[0:.01:12];
ex=i;
ind=1;
for x=i    
    if x<a1
        w1=-p1*x^2/(6*E*I)*(3*a1-x);
        M1=p1*a1*(1-x/a1);
    else
        w1=-p1*a1^2/(6*E*I)*(3*x-a1);
        M1=0;
    end
    
    if x<a2
        w2=-p2*x^2/(6*E*I)*(3*a2-x);
        M2=p2*a2*(1-x/a2);
    else
        w2=-p2*a2^2/(6*E*I)*(3*x-a2);
            M2=0;
    end
    
    w3=-m*x^2/(2*E*I);
    M3=m;    
    w4=-p4*x^2*(3*L-x)/(6*E*I);
    M4=p4*(L-x);
    
    if x<a3
        w5=-p3*x^2*(6*a3^2-4*x*a3+x^2)/(24*E*I);
        M5=1/2*p3*(a3-x)^2;
    else
        w5=   -p3*a3^3*(4*x-a3)/(24*E*I);
        M5=0;
    end
    if x<=12
        S1=p4;
    else
        S1=0;
    end
    if x<=8
        S2=p2;
        S3=(8-x)*p3;
    else
        S2=0;            
        S3=0;
    end
    if x<4
        S4=p1;
    else
        S4=0;
    end
    w(ind)=w1+w2+w3+w4+w5;
    M(ind)=M1+M2+M3+M4+M5;
    S(ind)=S1+S2+S3+S4;
    ind=ind+1;                        
end
M=-M;
S=-S;


% HW5 Q3 Exact ************************************************************
L = 1;
EI = 2e7;
p = 2000;
i = 0:0.1:L;

R2 = 5/8 * p * L;
R1 = 3/8 * p * L;

S = R1 - p .* (L-i);
S = -S;
M = R1 * (L-i) - 0.5 * p * (L-i).^2;
M = -M;
w = (p * (L-i)) / (48 * EI) .* (L^3 - 3 * L * (L-i).^2 + 2 * (L-i).^3);

ex = i;

% *************************************************************************

% HW6 Q2
L = 5; q = 5000; p = -4.9306e+03; b = 2*L/3;
E=2e5; I=1;
ex_1= 0:.01:L/3;
ex_2= L/3:0.01:L;
w_dist_a = 1/(E*I)*(q*ex_1.^5/(120*L) - q*L^3*ex_1/24 + q*L^4/30);
w_point_a_1 = p*b^2/(6*E*I)*(3*L - 3*ex_1 -b);
w_fin_1 = w_dist_a + w_point_a_1; 
w_dist_b = 1/(E*I)*(q*ex_2.^5/(120*L) - q*L^3*ex_2/24 + q*L^4/30);
w_dist_b_2 = p*(L-ex_2).^2/(6*E*I).*(3*b - L + ex_2);
w_fin_2 = w_dist_b + w_dist_b_2; 
w_analy = [w_fin_1 w_fin_2];
w = w_analy;
S_1 = -q*ex_1.^2/(2*L);
S_a_1 = -q*ex_2.^2/(2*L);
S_a_2 = -p;
S_a = S_a_1 + S_a_2;
S = [S_1  S_a];
M_1 = q*ex_1.^3/(6*L);
M_a_1 = q*ex_2.^3/(6*L);
M_a_2 = p*(ex_2 - L/3);
M_a = M_a_1 + M_a_2;
M = [M_1 M_a];
ex = [ex_1 ex_2];

