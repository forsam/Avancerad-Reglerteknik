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