%% For Question 4 part 3
K_small = [1 -1;
           -1 1];
re_1 = [1 0 0 0 0 0;
        0 0 0 1 0 0];
re_2 = [0 1 0 0 0 0;
        0 0 0 0 1 0];
re_3 = [sqrt(2)/2 sqrt(2)/2 0 0 0 0;
        0 0 0 sqrt(2)/2 sqrt(2)/2 0];
format short
disp(re_1' * K_small * re_1);
disp(re_2' * K_small * re_2);
disp(re_3' * K_small * re_3);

%% For solid mech
F = [sqrt(2), .75*sqrt(2), 0;
     -1, .75, sqrt(2)/4;
     1, -0.75, sqrt(2)/4];
C = F' * F;
disp(C);
x = [1; 1; 1];

disp(F*x);