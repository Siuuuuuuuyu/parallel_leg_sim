clc;
clear;

%% 定义腿部常量
L1a = 0.15/2;
L1u = 0.15;
L1d = 0.27;
m1u = 0.2;
m1d = 0.3;
L2a = 0.15/2;
L2u = 0.15;
L2d = 0.27;
m2u = 0.2;
m2d = 0.3;
lamda_1u = 0.5;
lamda_2u = 0.5;
lamda_1d = 0.5;
lamda_2d = 0.5;
I1u = 0.000375;
I2u = 0.000375;
I1d = 0.0018225;
I2d = 0.0018225;

%% 定义机体常量
ml = m1u + m2u + m1d + m2d; % 腿质量 kg
Rw = 0.06;                  % 驱动轮半径
Rl = 0.5 / 2;               % 驱动轮距离除以2
Lc = 0.025;                 % 机体质心到腿关节xz平面距离 m (这里机体质心与腿部关节连线与z轴平行）
mw = 0.668;                 % 驱动轮质量
mb = 12;                    % 机体质量 kg
Iw=0.0004656*3;             % 驱动轮转动惯量
% Ib=0.372204084;           % 机体俯仰转动惯量
Ib=0.3672204084;            % 机体俯仰转动惯量
Iz=0.328078099;             % 机体偏航转动惯量
g=9.81;                     % 重力加速度

%% 腿运动学逆解
syms l theta phi1 phi2;
l1 = sqrt(L1a^2 + l^2 + 2 * L1a * l * sin(theta));
l2 = sqrt(L2a^2 + l^2 - 2 * L2a * l * sin(theta));
phi1_ik = acos((L1u^2 + l1^2 - L1d^2) / (2 * L2u * l1)) ...
        +acos((l1^2 + L1a^2 - l^2) / (2 * L1a * l1));
phi2_ik = acos((L2u^2 + l2^2 - L2d^2) / (2 * L2u * l2)) ...
        +acos((l2^2 + L2a^2 - l^2) / (2 * L2a * l2));

%% 腿运动学正解
x_B1 = - L1u * cos(phi1);
y_B1 = - L1u * sin(phi1);
x_B2 = - (L1a + L2a) + L2u * cos(phi2);
y_B2 = - L2u * sin(phi2);
a_0 = x_B1 - x_B2;
b_0 = y_B2 - y_B1;
phi3 = asin((L2d^2 - a_0^2 - b_0^2 - L1d^2) / (sqrt((2 * b_0 * L1d)^2 + (2 * a_0 * L1d)^2))) ...
        + atan2(a_0, b_0);
x_C = x_B1 - L1d * cos(phi3);
y_C = y_B1 - L1d * sin(phi3);
l_phi = sqrt((x_C + L1a)^2 + y_C^2);
theta_phi = asin((- L1a - x_C) / l_phi);

%% 腿质心坐标解算
xT = - L1a - L2a;
yT = 0;
xB1 = - L1u * cos(phi1);
yB1 = - L1u * sin(phi1);
xB2 = - (L1a + L2a) + L2u * cos(phi2);
yB2 = - L2u * sin(phi2);
xC = - L1a - l_phi * sin(theta_phi);
yC = - l_phi * cos(theta_phi);
x1u = (xB1 - 0) * lamda_1u + 0;
y1u = (yB1 - 0) * lamda_1u + 0;
x1d = (xC - xB1) * lamda_1d + xB1;
y1d = (yC - yB1) * lamda_1d + yB1;
x2u = (xB2 - xT) * lamda_2u + xT;
y2u = (yB2 - yT) * lamda_2u + yT;
x2d = (xC - xB2) * lamda_2d + xB2;
y2d = (yC - yB2) * lamda_2d + yB2;
xmc = (x1u * m1u / ml + x1d * m1d / ml + x2u * m2u / ml + x2d * m2d / ml);
ymc = (y1u * m1u / ml + y1d * m1d / ml + y2u * m2u / ml + y2d * m2d / ml);

%% 腿部动力学参数解算
lamda_leg = ymc / y_C;
% 驱动轮到腿质心距离
Lw = (1 - lamda_leg) * l;
Lw = subs(Lw, [phi1, phi2], [phi1_ik, phi2_ik]);
Dmc_1u = sqrt((xmc - x1u)^2 + (ymc - y1u)^2);
Dmc_2u = sqrt((xmc - x2u)^2 + (ymc - y2u)^2);
Dmc_1d = sqrt((xmc - x1d)^2 + (ymc - y1d)^2);
Dmc_2d = sqrt((xmc - x2d)^2 + (ymc - y2d)^2);
% 腿转动惯量
Il = I1u + I2u + I1d + I2d + m1u * Dmc_1u^2 + m2u * Dmc_2u^2 + m1d * Dmc_1d^2 + m2d * Dmc_2d^2;
Il = subs(Il, [phi1, phi2], [phi1_ik, phi2_ik]);

%% 多项式拟合得到状态反馈矩阵K
LK_pair_array=struct();
x_ll = [];
y_lr = [];
z_k = [];
for i = 1:6
    Ll = 0.05 + 0.05 * i;
    for j = 1:6
        Lr = 0.05 + 0.05 * j;
% for i = 1
%     Ll = 0.11;
%     for j = 1
%         Lr = 0.11;
        Lwl = double(subs(Lw, [l, theta], [Ll, 0]));    % 左驱动轮到左腿质心距离
        Lbl = Ll - Lwl;                                 % 左腿关节到左腿质心距离
        Lwr = double(subs(Lw, [l, theta], [Lr, 0]));    % 右驱动轮到右腿质心距离
        Lbr = Lr - Lwr;                                 % 右腿关节到右腿质心距离
        Ill = double(subs(Il, [l, theta], [Ll, 0]));    % 左腿转动惯量
        Ilr = double(subs(Il, [l, theta], [Lr, 0]));    % 右腿转动惯量
        K = get_K_LQR(Ll, Lr, Lwl, Lbl, Lwr, Lbr, Lc, Ill, Ilr, ml, Rw, Rl, mw, mb, Iw, Ib, Iz);
        LK_pair_array((i - 1) * 6 + j).Ll = Ll;
        LK_pair_array((i - 1) * 6 + j).Lr = Lr;
        LK_pair_array((i - 1) * 6 + j).K = K;
        x_ll = [x_ll;Ll];
        y_lr = [y_lr;Lr];
        z_k = [z_k;K];
    end
end
n = size(x_ll);
x_poly = [ones(n), x_ll, y_lr, x_ll.*y_lr, x_ll.^2, y_lr.^2];
[row, column] = size(K);
% K = K1 + K2 * Ll + K3 * Lr + K4 * Ll * Lr + K5 * Ll^2 + K6 * Lr^2
Kk = struct();
K1 = zeros(row, column);
K2 = zeros(row, column);
K3 = zeros(row, column);
K4 = zeros(row, column);
K5 = zeros(row, column);
K6 = zeros(row, column);
Kk.Kk1 = K1;
Kk.Kk2 = K2;
Kk.Kk3 = K3;
Kk.Kk4 = K4;
Kk.Kk5 = K5;
Kk.Kk6 = K6;
for i = 1:row
    for j = 1:column
        z = [];
        for c = 1:n(1, 1)    
            z = [z;z_k(i + row * (c - 1), j)];
        end
        % 多项式回归
        b = regress(z, x_poly);
        Kk.Kk1(i, j) = b(1, 1);
        Kk.Kk2(i, j) = b(2, 1);
        Kk.Kk3(i, j) = b(3, 1);
        Kk.Kk4(i, j) = b(4, 1);
        Kk.Kk5(i, j) = b(5, 1);
        Kk.Kk6(i, j) = b(6, 1);
    end
end
