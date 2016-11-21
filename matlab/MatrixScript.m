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
     
% for PID
Beta = [0;b1;0;b2];
% This is for ss-simulation
BetaSim = [0, 0;b1, l_w;0, 0;b2, l_b];

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
P3 = -1;

% Gnum(3) = A
% Gnum(4) = B
% Gden(2) = C
% Gden(3) = D
% Gden(4) = E

Kd = (-(P1 + P2 + P3) - Gden(2))/Gnum(3);
Kp = (P1*P2 + P1*P3 + P2*P3 - Gden(3))/(Gnum(3));
Ki = (-P1*P2*P3 - Gden(4))/(Gnum(3));
D = tf([Kd Kp Ki],[1 0]);

T = feedback(D*G,1);


%%förfilter
%[Tnum,Tden]=tfdata(T);
%F = tf([Tden{1}(4)],Gnum(3)*[Kd,Kp,Ki]);
%T = minreal(F*T);



