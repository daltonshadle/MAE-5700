
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