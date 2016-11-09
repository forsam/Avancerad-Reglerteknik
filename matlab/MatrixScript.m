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

P1 = 100;
P2 = 200;
P3 = 300;
P4 = 800;

% Gden(2) = C
% Gden(3) = D
% Gden(4) = E
% Gnum(3) = A
% Gnum(4) = B

Kd = (P1 + P2 + P3 + P4 + Gden(2))/Gnum(3);
Kp = (P1*P2 + P1*P3 + P1*P4 + P2*P3 + P2*P4 + P3*P4 - Kd*Gnum(4) + Gden(3))/(Gnum(3));
Ki = (P1*P2*P3 + P1*P2*P4 + P1*P3*P4 + P2*P3*P4 - Kp*Gnum(4) - Gden(4))/(Gnum(3));

D = tf([Kd Kp Ki],[1 0]);

T = D*G/(D*G + 1);
step(T)

