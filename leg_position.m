function [l0, phi0] = leg_position(phi1,phi4)
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
    Xc = Xb + l2 * cos(phi2);
    Yc = Yb + l2 * sin(phi2);

    l0 = sqrt((Xc - l5/2)^2 + Yc^2);
    phi0 = atan2(Yc, (Xc - l5/2));
end
