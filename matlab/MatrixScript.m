CopiedCode
% make the matrices
Gamma = [1, 0, 0, 0;
         0, G11, 0, G12;
         0, 0, 1, 0;
         0, G21, 0, G22];

Alpha = [0, 1, 0, 0;
         a11,a12,a13,a14;
         0, 0, 0, 1;
         a21,a22,a23,a24];
     
Beta = [0;b1;0;b2];

A = Gamma\Alpha;
B = Gamma\Beta;
C = [0,0,1,0];

% Solve this, make the transfer function which is
% Y(s) = G(s)*U(s)
% G(s) = Y(s)/U(s)

G = tf(ss(A,B,C,0));
[Gnum,Gden] = ss2tf(A,B,C,0);

