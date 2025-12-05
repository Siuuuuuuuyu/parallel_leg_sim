%%定义符号变量
syms phi1 phi2 phi3 phi4;
syms d_phi1 d_phi4;
syms Xc Yc Xb Yb Xd Yd;
syms F Tp;
l1 = 0.15;
l2 = 0.25;
l3 = l2;
l4 = l1;
l5 = 0.108;

Xb = l1 * cos(phi1);
Yb = l1 * sin(phi1);
Xd = l4 * cos(phi4) + l5;
Yd = l4 * sin(phi4);

lbd = sqrt((Xd - Xb)^2 + (Yd - Yb)^2);
A0 = 2 * l2 * (Xd - Xb);
B0 = 2 * l2 * (Yd - Yb);
C0 = l2^2 + lbd^2 - l3^2;

phi2 = 2 * atan2((B0 + sqrt(A0^2 + B0^2 - C0^2)), (A0 + C0));
phi3 = atan2((Yb - Yd + l2 * sin(phi2)), (Xb - Xd + l2 * cos(phi2)));

Xc = Xb + l2 * cos(phi2);
Yc = Yb + l2 * sin(phi2);

l0 = sqrt((Xc - l5/2)^2 + Yc^2);
phi0 = atan2(Yc, (Xc - l5/2));
position = [l0; phi0];

%%l0 = f(phi1, phi4)
%%phi0 = g(phi, phi4)
JacobianMatrix = jacobian([l0; phi0], [phi1; phi4]);
%%直接通过diff运算得到的结果不正确
% J11 = diff(l0, phi1);
% J12 = diff(l0, phi4);
% J21 = diff(phi0, phi1);
% J22 = diff(phi0, phi4);
% JacobianMatrix = [J11 J12;
%                   J21 J22];

% R = [cos(phi0-pi/2) -sin(phi0-pi/2);
%      sin(phi0-pi/2)  cos(phi0-pi/2)];
% M = [0 -1/l0;
%      1     0];

j11 = (l1 * sin(phi0 - phi3) * sin(phi1 - phi2)) / sin(phi3 - phi2);
j12 = (l1 * cos(phi0 - phi3) * sin(phi1 - phi2)) / (l0 * sin(phi3 - phi2));
j21 = (l4 * sin(phi0 - phi2) * sin(phi3 - phi4)) / sin(phi3 - phi2);
j22 = (l4 * cos(phi0 - phi2) * sin(phi3 - phi4)) / (l0 * sin(phi3 - phi2));
J = [j11 j12;
     j21 j22];

% torque = J * [F; Tp];
% matlabFunction(torque, "File", "leg_torque");
