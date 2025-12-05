%%腿长L0变化范围
L0s = 0.08:0.01:0.40;

%%根据L0计算出对应lqr增益矩阵
K_lqr = zeros(2, 6, length(L0s));

for i = 1:length(L0s)
    %%定义符号变量
    syms theta d_theta dd_theta;
    syms phi d_phi dd_phi;
    syms x d_x dd_x;
    syms T Tp N P Nf Nm Pm;

    %%结构参数
    mw = 0.5;                       %驱动轮质量
    R = 0.08;                       %驱动轮半径
    l = 0.02;                       %机体质心到转轴距离
    mp = 1;                         %摆杆质量
    M = 10;                         %机体质量
    g = 9.8;
    L = L0s(i) / 2;                 %摆杆质心到驱动轮距离
    Lm = L0s(i) / 2;                %摆杆质心到机体转轴距离
    Ip = (mp * (L+Lm)^2) / 12;      %摆杆绕质心转动惯量
    Iw = (mw * R^2) / 2;            %驱动轮转动惯量
    Im = (M * 0.1625) / 12;         %机体绕质心转动惯量 
    
    xb = x + (L+Lm) * sin(theta);
    d_xb = d_x + (L+Lm) * d_theta * cos(theta);
    dd_xb = dd_x + (L+Lm) * dd_theta * cos(theta) - (L+Lm) * d_theta^2 *sin(theta);

    Nm = M * (dd_xb - l * dd_phi * cos(phi) + l * d_phi^2 * sin(phi));
    Pm = M * (g - (L+Lm) * dd_theta * sin(theta) - (L+Lm) * d_theta^2 * cos(theta) - l * dd_phi * sin(phi) - l * d_phi^2 * cos(phi));

    N = Nm + mp * (dd_x + L * dd_theta * cos(theta) - L * d_theta^2 * sin(theta));
    P = Pm + mp * (g - L * d_theta^2 * cos(theta) - L * dd_theta * sin(theta));

    equ1 = dd_x - ((T - N * R) / ((Iw / R) + (mw * R)));
    equ2 = ((P * L) + (Pm * Lm)) * sin(theta) - ((N * L) + (Nm * Lm)) * cos(theta) + Tp - T - Ip * dd_theta;
    equ3 = Nm * l * cos(phi) + Pm * l * sin(phi) + Tp - Im * dd_phi;

    [dd_x, dd_theta, dd_phi] = solve(equ1, equ2, equ3, dd_x, dd_theta, dd_phi);

    JA = jacobian([d_theta, dd_theta, d_x, dd_x, d_phi, dd_phi], [theta, d_theta, x, d_x, phi, d_phi]);
    JB = jacobian([d_theta, dd_theta, d_x, dd_x, d_phi, dd_phi], [T, Tp]);
    
    A = vpa(subs(JA, [theta; d_theta; xb; d_xb; phi; d_phi], [0; 0; 0; 0; 0; 0]));
    B = vpa(subs(JB, [theta; d_theta; xb; d_xb; phi; d_phi], [0; 0; 0; 0; 0; 0]));

    %%转换为双精度数值矩阵
    A_num = double(A);
    B_num = double(B);

    %%判断能控性
    C = ctrb(A_num, B_num);
    r = rank(C);

    %%离散化处理
    [a, b] = c2d(eval(A), eval(B), 0.001);
    
    %%theta, d_theta, x, d_x, phi, d_phi 
    Q = diag([1 1 10 100 300 1]);
    %%T; Tp
    R = diag([1 0.5]);

    K_lqr(:, :, i) = dlqr(a, b, Q, R);
end

%%对K_lqr中每个元素进行拟合
K = sym("K", [2 6]);
syms L0;

for x = 1:2
    for y = 1:6
        p = polyfit(L0s, reshape(K_lqr(x, y, :), 1, length(L0s)), 3);
        K(x, y) = p(1) * L0^3 + p(2) * L0^2 + p(3) * L0 + p(4);
    end
end

matlabFunction(K, "File", "K_lqr");
