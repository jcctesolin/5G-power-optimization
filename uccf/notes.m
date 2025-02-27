
B = [1 2;
     3 4;
     5 6;
     7 8];

C = [1 2;
     0 4;
     5 0;
     0 0];

A = C * pinv(B);

disp(A);


A_line = A;

disp(A_line*B);