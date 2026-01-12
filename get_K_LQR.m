function K = get_K_LQR(Ll, Lr, Lwl, Lbl, Lwr, Lbr, Lc, Ill, Ilr, ml, Rw, Rl, mw, mb, Iw, Ib, Iz)
% Lwl, Lwr : 驱动轮到腿部质心距离
% Lbl, Lbr : 机体到腿部质心距离
% Lc : 机体质心到髋关节距离
% Rw : 驱动轮半径
% Rl : 驱动轮轮距
syms s d_s dd_s;                        % 自然坐标系下机器人水平方向移动距离
syms phi d_phi dd_phi;                  % 机体偏航角yaw
syms theta_ll d_theta_ll dd_theta_ll;
syms theta_lr d_theta_lr dd_theta_lr;   % 腿倾角theta_l
syms theta_b d_theta_b dd_theta_b;      % 机体倾角
syms Twl Twr                            % 驱动轮力矩
syms Tbl Tbr;                           % 髋关节输出力矩
syms theta_wl d_theta_wl dd_theta_wl;
syms theta_wr d_theta_wr dd_theta_wr;   % 驱动轮转角theta_w
syms sb d_sb dd_sb;                     % 自然坐标系下机体质心水平移动距离
syms hb d_hb dd_hb;                     % 自然坐标系下机体质心竖直移动距离
syms sll d_sll dd_sll slr d_slr dd_slr; % 自然坐标系下左右腿质心水平移动距离
syms hll d_hll dd_hll hlr d_hlr dd_hlr; % 自然坐标系下左右腿质心竖直移动距离
syms fl  fr;                            % 左右轮所受摩擦力
syms Fwsl Fwsr;                         % 驱动轮对腿水平方向作用力
syms Fwhl Fwhr;                         % 驱动轮对腿竖直方向作用力
syms Fbsl Fbsr;                         % 腿对机体水平方向的作用力
syms Fbhl Fbhr;                         % 腿对机体竖直方向的作用力
syms d2_theta_wl d2_theta_wr d2_theta_ll d2_theta_lr d2_theta_b;
g = 9.81;

% 运动学
% s = (Rw * (theta_wl + theta_wr)) / 2;
% sll = Rw * theta_wl + Lwl * sin(theta_ll);
% slr = Rw * theta_wr + Lwr * sin(theta_lr);
% hb = (Ll * cos(theta_ll) + Lr * cos(theta_lr)) / 2 + Lc * cos(theta_b);
% hll = hb - Lbl * cos(theta_ll) - Lc * cos(theta_b);
% hlr = hb - Lbr * cos(theta_lr) - Lc * cos(theta_b);

% Rw * theta_wl = (sb - Lc * sin(theta_b)) - Rl * phi - Ll * sin(theta_ll);
% Rw * theta_wr = (sb - Lc * sin(theta_b)) + Rl * phi - Lr * sin(theta_lr);
% 联立上式得
% phi = (Rw * (-theta_wl + theta_wr) - Ll * sin(theta_ll) + Lr * sin(theta_lr)) / (2 * Rl);
% sb = (Rw * (theta_wl + theta_wr) + Ll * sin(theta_ll) + Lr * sin(theta_lr)) / 2 + Lc * sin(theta_b);

% d_s = Rw / 2 * (d_theta_wl + d_theta_wr);
% dd_s = Rw / 2 * (dd_theta_wl + dd_theta_wr);
% d_phi = (Rw * (-d_theta_wl + d_theta_wr) ...
%         - Ll * cos(theta_ll) * d_theta_ll ...
%         + Lr * cos(theta_lr) * d_theta_lr) / (2 * Rl);
dd_phi = (Rw * (-dd_theta_wl + dd_theta_wr) ...
        - Ll * cos(theta_ll) * dd_theta_ll ...
        + Ll * sin(theta_ll) * (d_theta_ll^2) ...
        + Lr * cos(theta_lr) * dd_theta_lr ...
        - Lr * sin(theta_lr) * (d_theta_lr^2)) / (2 * Rl);
% d_sb = (Rw * (d_theta_wl + d_theta_wr) ...
%         + Ll * cos(theta_ll) * d_theta_ll ...
%         + Lr * cos(theta_lr) * d_theta_lr) / 2 ...
%         + Lc * cos(theta_b) * d_theta_b;
dd_sb = (Rw * (dd_theta_wl + dd_theta_wr) ...
        + Ll * cos(theta_ll) * dd_theta_ll ...
        - Ll * sin(theta_ll) * (d_theta_ll^2) ...
        + Lr * cos(theta_lr) * dd_theta_lr ...
        - Lr * sin(theta_lr) * (d_theta_lr^2)) / 2 ...
        + Lc * cos(theta_b) * dd_theta_b ...
        - Lc * sin(theta_b) * (d_theta_b^2);
% d_hb = (- Ll * sin(theta_ll) * d_theta_ll ...
%         - Lr * sin(theta_lr) * d_theta_lr) / 2 ...
%         - Lc * sin(theta_b) * d_theta_b;
dd_hb = (- Ll * sin(theta_ll) * dd_theta_ll ...
        - Ll * cos(theta_ll) * (d_theta_ll^2) ...
        - Lr * sin(theta_lr) * dd_theta_lr ...
        - Lr * cos(theta_lr) * (d_theta_lr^2)) / 2 ...
        - Lc * sin(theta_b) * dd_theta_b ...
        - Lc * cos(theta_b) * (d_theta_b^2);

% d_sll = Rw * d_theta_wl + Lwl * cos(theta_ll) * d_theta_ll;
dd_sll = Rw * dd_theta_wl + Lwl * (cos(theta_ll) * dd_theta_ll - sin(theta_ll) * (d_theta_ll^2));
% d_slr = Rw * d_theta_wr + Lwr * cos(theta_lr) * d_theta_lr;
dd_slr = Rw * dd_theta_wr + Lwr * (cos(theta_lr) * dd_theta_lr - sin(theta_lr) * (d_theta_lr^2));

% d_hll = d_hb + Lbl * sin(theta_ll) * d_theta_ll + Lc * sin(theta_b) * dd_theta_b;
dd_hll = dd_hb + Lbl * (sin(theta_ll) * dd_theta_ll ...
        + cos(theta_ll) * (d_theta_ll^2)) ...
        + Lc * sin(theta_b) * dd_theta_b ...
        + Lc * cos(theta_b) * (d_theta_b^2);
% d_hlr = d_hb + Lbr * sin(theta_lr) * d_theta_lr + Lc * sin(theta_b) * dd_theta_b;
dd_hlr = dd_hb + Lbr * (sin(theta_lr) * dd_theta_lr ...
        + cos(theta_lr) * (d_theta_lr^2)) ...
        + Lc * sin(theta_b) * dd_theta_b ...
        + Lc * cos(theta_b) * (d_theta_b^2);

% 动力学-牛顿欧拉法
% 对驱动轮
% mw * Rw * dd_theta_wl = fl - Fwsl;
% Iw * dd_theta_wl = Twl - fl * Rw;
% mw * Rw * dd_theta_wr = fr - Fwsr;
% Iw * dd_theta_wr = Twr - fr * Rw;
% 对腿
% ml * dd_sll = Fwsl - Fbsl;
% ml * dd_hll = Fwhl - Fbhl - ml *g;
% Ill * dd_theta_ll = (Fwhl * Lwl + Fbhl * Lbl) * sin(theta_ll) - (Fwsl * Lwl + Fbsl * Lbl) * cos(theta_ll) - Twl + Tbl;
% ml * dd_slr = Fwsr - Fbsr;
% ml * dd_hlr = Fwhr - Fbhr - ml *g;
% Ilr * dd_theta_lr = (Fwhr * Lwr + Fbhr * Lbr) * sin(theta_lr) - (Fwsr * Lwr + Fbsr * Lbr) * cos(theta_lr) - Twr + Tbr;
% 对机体
% mb * dd_sb = Fbsl + Fbsr;
% mb * dd_hb = Fbhl + Fbhr - mb *g;
% Ib * dd_theta_b = (Fbhl + Fbhr) * Lc *sin(theta_b) - (Fbsl + Fbsr) * Lc * cos(theta_b) - (Tbl + Tbr);
% 对机体偏航角
% Iz * dd_phi = (fr - fl) * Rl;
% 左右腿支持力相等
% Fwhl = Fwhr;

% 求解约束力
Fbhl = (mb * dd_hb + mb * g - ml * dd_hll + ml * dd_hlr) / 2;
Fbhr = (mb * dd_hb + mb * g + ml * dd_hll - ml * dd_hlr) / 2;
Fwhl = ml * dd_hll + ml * g + Fbhl;
Fwhr = ml * dd_hlr + ml * g + Fbhr;
Fbsl = (mb * dd_sb - Iz / Rl * dd_phi + mw * Rw * (dd_theta_wr - dd_theta_wl) + ml * (dd_slr - dd_sll)) / 2;
Fbsr = (mb * dd_sb + Iz / Rl * dd_phi - mw * Rw * (dd_theta_wr - dd_theta_wl) - ml * (dd_slr - dd_sll)) / 2;
Fwsl = ml * dd_sll + Fbsl;
Fwsr = ml * dd_slr + Fbsr;
fl = mw * Rw * dd_theta_wl + Fwsl;
fr = mw * Rw * dd_theta_wr + Fwsr;

% 未使用的动力学方程
equ1 = Iw * dd_theta_wl == Twl - fl * Rw;
equ2 = Iw * dd_theta_wr == Twr - fr * Rw;
equ3 = Ill * dd_theta_ll == (Fwhl * Lwl + Fbhl * Lbl) * sin(theta_ll) ...
                            - (Fwsl * Lwl + Fbsl * Lbl) * cos(theta_ll) ...
                            - Twl + Tbl;
equ4 = Ilr * dd_theta_lr == (Fwhr * Lwr + Fbhr * Lbr) * sin(theta_lr) ...
                            - (Fwsr * Lwr + Fbsr * Lbr) * cos(theta_lr) ...
                            - Twr + Tbr;
equ5 = Ib * dd_theta_b == (Fbhl + Fbhr) * Lc *sin(theta_b) ...
                            - (Fbsl + Fbsr) * Lc * cos(theta_b) ...
                            - (Tbl + Tbr);
dx_part = [dd_theta_wl dd_theta_wr dd_theta_ll dd_theta_lr dd_theta_b];
fx = [equ1 equ2 equ3 equ4 equ5];
[d2_theta_wl, d2_theta_wr, d2_theta_ll, d2_theta_lr, d2_theta_b] = solve(fx, dx_part);

% x = [s; d_s; phi; d_phi; theta_ll; d_theta_ll; theta_lr; d_theta_lr; theta_b; d_theta_b]
% u = [Twl; Twr; Tbl; Tbr]

% A矩阵中的特定项
a25 = (Rw / 2) * (diff(d2_theta_wl, theta_ll) + diff(d2_theta_wr, theta_ll));
a27 = (Rw / 2) * (diff(d2_theta_wl, theta_lr) + diff(d2_theta_wr, theta_lr));
a29 = (Rw / 2) * (diff(d2_theta_wl, theta_b) + diff(d2_theta_wr, theta_b));
a45 = (Rw / (2 * Rl)) * (- diff(d2_theta_wl, theta_ll) + diff(d2_theta_wr, theta_ll)) ...
        - (Ll / (2 * Rl)) * diff(d2_theta_ll, theta_ll) ...
        + (Lr / (2 * Rl)) * diff(d2_theta_lr, theta_ll);
a47 = (Rw / (2 * Rl)) * (- diff(d2_theta_wl, theta_lr) + diff(d2_theta_wr, theta_lr)) ...
        - (Ll / (2 * Rl)) * diff(d2_theta_ll, theta_lr) ...
        + (Lr / (2 * Rl)) * diff(d2_theta_lr, theta_lr);
a49 = (Rw / (2 * Rl)) * (- diff(d2_theta_wl, theta_b) + diff(d2_theta_wr, theta_b)) ...
        - (Ll / (2 * Rl)) * diff(d2_theta_ll, theta_b) ...
        + (Lr / (2 * Rl)) * diff(d2_theta_lr, theta_b);
a65 = diff(d2_theta_ll, theta_ll);
a67 = diff(d2_theta_ll, theta_lr);
a69 = diff(d2_theta_ll, theta_b);
a85 = diff(d2_theta_lr, theta_ll);
a87 = diff(d2_theta_lr, theta_lr);
a89 = diff(d2_theta_lr, theta_b);
a105 = diff(d2_theta_b, theta_ll);
a107 = diff(d2_theta_b, theta_lr);
a109 = diff(d2_theta_b, theta_b);

% B矩阵中的特定项
b21 = (Rw / 2) * (diff(d2_theta_wl, Twl) + diff(d2_theta_wr, Twl));
b22 = (Rw / 2) * (diff(d2_theta_wl, Twr) + diff(d2_theta_wr, Twr));
b23 = (Rw / 2) * (diff(d2_theta_wl, Tbl) + diff(d2_theta_wr, Tbl));
b24 = (Rw / 2) * (diff(d2_theta_wl, Tbr) + diff(d2_theta_wr, Tbr));
b41 = (Rw / (2 * Rl)) * (- diff(d2_theta_wl, Twl) + diff(d2_theta_wr, Twl)) ...
        - (Ll / (2 * Rl)) * diff(d2_theta_ll, Twl) ...
        + (Lr / (2 * Rl)) * diff(d2_theta_lr, Twl);
b42 = (Rw / (2 * Rl)) * (- diff(d2_theta_wl, Twr) + diff(d2_theta_wr, Twr)) ...
        - (Ll / (2 * Rl)) * diff(d2_theta_ll, Twr) ...
        + (Lr / (2 * Rl)) * diff(d2_theta_lr, Twr);
b43 = (Rw / (2 * Rl)) * (- diff(d2_theta_wl, Tbl) + diff(d2_theta_wr, Tbl)) ...
        - (Ll / (2 * Rl)) * diff(d2_theta_ll, Tbl) ...
        + (Lr / (2 * Rl)) * diff(d2_theta_lr, Tbl);
b44 = (Rw / (2 * Rl)) * (- diff(d2_theta_wl, Tbr)+diff(d2_theta_wr, Tbr)) ...
        - (Ll / (2 * Rl)) * diff(d2_theta_ll, Tbr) ...
        + (Lr / (2 * Rl)) * diff(d2_theta_lr, Tbr);
b61 = diff(d2_theta_ll, Twl);
b62 = diff(d2_theta_ll, Twr);
b63 = diff(d2_theta_ll, Tbl);
b64 = diff(d2_theta_ll, Tbr);
b81 = diff(d2_theta_lr, Twl);
b82 = diff(d2_theta_lr, Twr);
b83 = diff(d2_theta_lr, Tbl);
b84 = diff(d2_theta_lr, Tbr);
b101 = diff(d2_theta_b, Twl);
b102 = diff(d2_theta_b, Twr);
b103 = diff(d2_theta_b, Tbl);
b104 = diff(d2_theta_b, Tbr);

A_sym = sym(zeros(10, 10));
B_sym = sym(zeros(10, 4));
A_sym(2, 5) = a25; A_sym(2, 7) = a27; A_sym(2, 9) = a29;
A_sym(4, 5) = a45; A_sym(4, 7) = a47; A_sym(4, 9) = a49;
A_sym(6, 5) = a65; A_sym(6, 7) = a67; A_sym(6, 9) = a69;
A_sym(8, 5) = a85; A_sym(8, 7) = a87; A_sym(8, 9) = a89;
A_sym(10, 5) = a105; A_sym(10, 7) = a107; A_sym(10, 9) = a109;

B_sym(2, 1) = b21; B_sym(2, 2) = b22; B_sym(2, 3) = b23; B_sym(2, 4) = b24;
B_sym(4, 1) = b41; B_sym(4, 2) = b42; B_sym(4, 3) = b43; B_sym(4, 4) = b44;
B_sym(6, 1) = b61; B_sym(6, 2) = b62; B_sym(6, 3) = b63; B_sym(6, 4) = b64;
B_sym(8, 1) = b81; B_sym(8, 2) = b82; B_sym(8, 3) = b83; B_sym(8, 4) = b84;
B_sym(10, 1) = b101; B_sym(10, 2) = b102; B_sym(10, 3) = b103; B_sym(10, 4) = b104;

A = zeros(10, 10);
B = zeros(10, 4);
x = [s d_s phi d_phi theta_ll d_theta_ll theta_lr d_theta_lr theta_b d_theta_b];
u = [Twl Twr Tbl Tbr];
x0 = [0 0 0 0 0 0 0 0 0 0];
u0 = [0 0 0 0];
A_sym = subs(A_sym, x, x0);
A_sym = subs(A_sym, u, u0);
A = double(A_sym);
A(1,2) = 1; A(3, 4) = 1; A(5, 6) = 1; A(7, 8) = 1; A(9, 10) = 1;

B_sym = subs(B_sym, x, x0);
B_sym = subs(B_sym, u, u0);
B = double(B_sym);

% x = [s; d_s; phi; d_phi; theta_ll; d_theta_ll; theta_lr; d_theta_lr; theta_b; d_theta_b]
Q = diag([500 300 300 30 600 80 600 80 30000 1000]);
% u = [Twl; Twr; Tbl; Tbr]
R = diag([25 25 3 3]);

K = lqr(A, B, Q, R);
end
