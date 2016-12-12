CopiedCode
 
% make the matrices
Gamma = [G11, 0, 0, 0;
         0, G22, 0, G24;
         0, 0, G33, 0;
         0, G42, 0, G44];
 
Alpha = [a11, a12, a13, a14;
         a21,a22,a23,a24;
         a31, a32, a33, a34;
         a41,a42,a43,a44];
     
% for PID
Beta = [0;b2;0;b4];
% This is for ss-simulation
BetaSim = [0, 0;b2, l_w;0, 0;b4, l_b];
 
% For calculations for the PID
A = Gamma\Alpha;
B = Gamma\Beta;
C = [0,0,1,0];
D = 0;
 
% For ss-simulation
Csim = eye(4);
Dsim = zeros(4,2);
Bsim = Gamma\BetaSim;
 
 
%make the transfer function
G = tf(ss(A,B,C,0));
[Gnum,Gden] = ss2tf(A,B,C,D);
%Cleaning the function!
G = tf([Gnum(3) 0],[1 Gden(2) Gden(3) Gden(4)]);
 
%pole placement
[p z] = pzmap(G);
P1 = p(1);
P2 = p(3);
P3 = -5;
 
% Gnum(3) = A
% Gnum(4) = B
% Gden(2) = C
% Gden(3) = D
% Gden(4) = E
 
Kd = (-(P1 + P2 + P3) - Gden(2))/Gnum(3);
Kp = (P1*P2 + P1*P3 + P2*P3 - Gden(3))/(Gnum(3));
Ki = (-P1*P2*P3 - Gden(4))/(Gnum(3));
D = tf([Kd Kp Ki],[1 0]);
 
% The Transfer function!
T = feedback(D*G,1);
T = minreal(T);
% The Bandwidth!
BW = bandwidth(D*G);
 
% Defined from the book! 
samplingTime = 1/(25*BW);
%TS = 1/(4*pi*bandwidth(T))
%%förfilter
%[Tnum,Tden]=tfdata(T);
%F = tf([Tden{1}(4)],Gnum(3)*[Kd,Kp,Ki]);
%T = minreal(F*T);
 
% Define controllability and oobservability
Controllability = ctrb(A,B);
Observability = obsv(A,C);
 
% test to create the New pole placements usins pole allocation with
% Controller!!
% Do the upper companion stuff
[Ac,Bc,Cc,Dc] = tf2ss(Gnum,Gden);
tn = [0 0 0 1]/Controllability;
Tinv = [tn*A^3;tn*A^2;tn*A;tn*eye(4)];
 
% bestämmer egenskaper för dominant andragrads system!!
overshoot = 0.1; % andel overshoot!!
ksi = sqrt((log(overshoot)^2)/(pi^2 + log(overshoot)^2)); % underdämpat
settlingTime = 2;
omega_n = 3.9/(ksi*settlingTime);
s1 = -ksi*omega_n + omega_n*sqrt(ksi^2 - 1);
s2 = -ksi*omega_n - omega_n*sqrt(ksi^2 - 1);
pol= -5;
pc = eig(A);
pc = [s1,s2,pc(2),pc(3)];
K = acker(A,B,pc);
 
%L = (acker(A',C',somePoles))';
 
%% LQR methods!!
C1 = [1 0 0 0;0 1 0 0; 0 0 sqrt(3) 0;0 0 0 1];
R = 1;
rho = 100;
Q = rho*diag(C1)*diag(C1)';
lqrK = lqr(A,B,Q,R);
 
% mimo system use place instead of acker!!
[G1num,G1den] = ss2tf(A,B,diag(C1)',0);
G1 = tf(G1num,G1den);
G2 = tf([G1num(1),-G1num(2),G1num(3),-G1num(4),G1num(5)],[G1den(1),-G1den(2),G1den(3),-G1den(4),G1den(5)]);
 
sysGG = G1*G2;
LQRRoots = rlocus(sysGG, 365);
stableRoots = LQRRoots(LQRRoots < 0);
 
K = place(A,B,stableRoots');
 
 
%% full order Luenberger observer
C = [1 0 0 0; 0 0 1 0];
L = (place(A',C',stableRoots'*6))';
 
%% reduced order Luenberger observer
Cacc = [1 0 0 0];
Cacc_ = [0 0 1 0];
V = [0 1 0 0;0 0 1 0;0 0 0 1];
invT = [Cacc;V];
T = inv(invT); 
 
Cy = Cacc_(1);
Cx = Cacc_(2:4);
 
By = B(1);
Bx = B(2:4);
 
Ayy = A(1,1);
Axx = A(2:4,2:4);
Axy = A(2:4,1);
Ayx = A(1,2:4);
 
Li = (place(Axx',[Ayx;Cx]',stableRoots(1:3)'*3))';
 
Lacc = Li(:,1);
Lacc_ = Li(:,2);
 
 
M1 = Axx - Lacc*Ayx - Lacc_*Cx;
M2 = Bx - Lacc*By;
M3 = Axy - Lacc*Ayy - Lacc_*Cy;
M4 = Lacc_;
M5 = Lacc;
M6 = T(1,:)';
M7 = T(1:4,2:4);
 
%discreteize system
 
discSys = c2d(ss(A,B,C,0),samplingTime);
Ad = discSys.A;
Bd = discSys.B;
Cd = discSys.C;
Dd = discSys.D;
 
z = exp(stableRoots*3*samplingTime);
 
Cdy = Cacc_(1);
Cdx = Cacc_(2:4);
 
Bdy = Bd(1);
Bdx = Bd(2:4);
 
Adyy = Ad(1,1);
Adxx = Ad(2:4,2:4);
Adxy = Ad(2:4,1);
Adyx = Ad(1,2:4);
 
Ldi = (place(Adxx',[Adyx;Cdx]',z(1:3)))';
 
Ldacc = Ldi(:,1);
Ldacc_ = Ldi(:,2);
 
 
Md1 = Adxx - Ldacc*Adyx - Ldacc_*Cdx;
Md2 = Bdx - Ldacc*Bdy;
Md3 = Adxy - Ldacc*Adyy - Ldacc_*Cdy;
Md4 = Ldacc_;
Md5 = Ldacc;
Md6 = T(1,:)';
Md7 = T(1:4,2:4);
 
kP = Kp;
kI = Ki;
kD = Kd;
D = 0;
 
Ld = (place(Ad',Cd',z'))';
Kd = place(Ad,Bd,z);
