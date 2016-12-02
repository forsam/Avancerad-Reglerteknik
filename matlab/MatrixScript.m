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

K = acker(A,B,pc)

%L = (acker(A',C',somePoles))';
C1 = [1 0 1 0];
R = 1;
rho = 100;
Q = rho*C1'*C1;
lqrK = lqr(A,B,Q,R)


% mimo system use place instead of acker!

syms s
N = det([eye(4)*s-Ac, -Bc ; C1,0]);
D = det([eye(4)*s-Ac]);
minusN = det(-[eye(4)*s-Ac, -Bc ; C1,0]);
minusD = det([-eye(4)*s-Ac]);

findK = minusD*D + rho*minusN*N;
koeff = sym2poly(findK);
theRoots = roots(koeff);
theStableRoots = theRoots( real(theRoots) < -0.0001 );


K = place(Ac,Bc,theStableRoots')
