function [K1, K2, K3, K4, K5] = getKBarRelations(r1, r2, r3, r4)
    K1 = r1/r2;
    K2 = r1/r4;
    K3 = (r2^2-r3^2+r4^2+r1^2)/(2*r2*r4);
    K4 = r1/r3;
    K5 = (r4^2-r1^2-r2^2-r3^2)/(2*r2*r3);
end
