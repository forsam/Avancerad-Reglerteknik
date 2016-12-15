clear all
close all
clc
g = 9.8;
b_f = 0;
R_m = 4.4;
L_m = 0;
b_m = 0;
K_e = 0.444;
K_t = 0.470;
m_b = 0.381;
l_b = 0.112;
I_b = 0.00616;
m_w = 0.036;
l_w = 0.021;
I_w = 0.00000746;
% Helping constant
C_1 = -(b_f + K_t*K_e/R_m);
% Fitting constants for matrices
% Gamma
G11 = 1;
G22 = I_w/l_w + m_w*l_w + m_b*l_w;
G24 = m_b*l_b*l_w;
G33 = 1;
G42 = m_b*l_b;
G44 = I_b + m_b*l_b^2;
% Alpha
a11 = 0;
a12 = 1;
a13 = 0;
a14 = 0;
a21 = 0;
a22 = C_1/l_w;
a23 = 0;
a24 = -C_1;
a31 = 0;
a32 = 0;
a33 = 0;
a34 = 1;
a41 = 0;
a42 = -C_1/l_w;
a43 = m_b*l_b*g;
a44 = C_1;
% Beta
b2 = K_t/R_m;
b4 = -b2;
% make the matrices
Gamma = [G11, 0, 0, 0;
         0, G22, 0, G24;
         0, 0, G33, 0;
         0, G42, 0, G44];

Alpha = [a11, a12, a13, a14;
         a21,a22,a23,a24;
         a31, a32, a33, a34;
         a41,a42,a43,a44];  
 
Beta = [0;b2;0;b4];

% The system
A = Gamma\Alpha;
B = Gamma\Beta;
C = [1,0,0,0; 0 0 1 0];
D = 0;

%% Discreteized system

Frequency = 100;
samplingPeriod = 1/Frequency;
discSys = c2d(ss(A,B,C,0),samplingPeriod);

Ad = discSys.A;
Bd = discSys.B;
Cd = discSys.C;
Dd = discSys.D;

%% Discrete LQR method

Qd_value = 50;
Qd = Qd_value*[12 0 0 0;0 3 0 0; 0 0 2 0;0 0 0 2];
Rd = 1;

[Kd,S,e] = dlqr(Ad,Bd,Qd,Rd);
chi = 6;

%% Discrete reduced order Luenberger observer
Cacc = C(1,:);
Cacc_ = C(2,:);
V = [0 1 0 0;0 0 1 0;0 0 0 1];
invT = [Cacc;V];
T = inv(invT); 

% Faster observer poles.
z = exp(chi*log(e));

Cdy = Cacc_(1);
Cdx = Cacc_(2:4);

Bdy = Bd(1);
Bdx = Bd(2:4);

Adyy = Ad(1,1);
Adxx = Ad(2:4,2:4);
Adxy = Ad(2:4,1);
Adyx = Ad(1,2:4);

Ldi = (place(Adxx',[Adyx;Cdx]',z(2:4)))';

Ldacc = Ldi(:,1);
Ldacc_ = Ldi(:,2);

Md1 = Adxx - Ldacc*Adyx - Ldacc_*Cdx;
Md2 = Bdx - Ldacc*Bdy;
Md3 = Adxy - Ldacc*Adyy - Ldacc_*Cdy;
Md4 = Ldacc_;
Md5 = Ldacc;
Md6 = T(1,:)';
Md7 = T(1:4,2:4);

% Full order discrete Luenberger observer
Ld = (place(Ad',Cd',z'))';

%% Calculate Nx and Nu
CompensatorMatrix = [(Ad - eye(4)), Bd; Cd(1,:), Dd(1)];
baseN = CompensatorMatrix\[0 0 0 0 1]'; 
Nxd = baseN(1:4);
Nud = baseN(5);

% The samplingperiod.
fSamplingPeriod = samplingPeriod;

% Do not modify these variables
iNumberOfEncoderSteps	= 720;
fGyroConversionFactor	= -1/131;
fWheelRadius			= 0.0216; % [m]
load('GyroBias.mat');
fGyroBias = -215;
