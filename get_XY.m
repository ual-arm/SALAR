function [XY] = get_XY(x, context)
    np = context.lenRefXY;
    r1=x(1); r2=x(2); r3=x(3); r4=x(4); rcx=x(5); rcy=x(6);

    switch context.prescribed
        case 1
            th2i = context.th2_n;
            if length(x)>6+1
                x0=x(6+1); y0=x(6+2); th0=x(6+3);
            else
                x0=0; y0=0; th0=0;
            end
        case 2
            th2i = x(7) + context.th2_n;
            if length(x)>6+1
                x0=x(7+1); y0=x(7+2); th0=x(7+3);
            else
                x0=0; y0=0; th0=0;
            end
        case 3
            th2i=x(7:6+np);
            if length(x)>6+np
                x0=x(6+np+1); y0=x(6+np+2); th0=x(6+np+3);
            else
                x0=0; y0=0; th0=0;
            end
    end
    
    XY=zeros(np,2);

    [K1, ~, ~, K4, K5] = getKBarRelations(r1, r2, r3, r4);
    for i=1:np
        th2=th2i(i);
        th3 = getTheta3(K1, K4, K5, th2);
        XY(i,:) = traceCXY(r2, rcx, rcy, th2, th3, x0, y0, th0);
    end
end
