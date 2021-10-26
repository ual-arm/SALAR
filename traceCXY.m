function pXY = traceCXY(r2, rcx, rcy, th2, th3, x0, y0, th0)
    Cxr = r2*cos(th2) + rcx*cos(th3) - rcy*sin(th3);
    Cyr = r2*sin(th2) + rcx*sin(th3) + rcy*cos(th3);
    
    CXY= [cos(th0) -sin(th0); sin(th0) cos(th0)]*[Cxr; Cyr]+[x0; y0];    
    pXY = [CXY(1),CXY(2)];
end
