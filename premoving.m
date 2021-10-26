function TPs = premoving(omega2, x, th2simul, calcs)
    alpha2 = 0;
    r1 = x(1); r2 = x(2); r3 = x(3); r4 = x(4); rcx = x(5); rcy =x(6);
    x0=x(7); y0=x(8); th0=x(9);
    if isnan(x0); x0=0; end
    if isnan(y0); y0=0; end
    if isnan(th0); th0=0; end

    xA=zeros(length(th2simul),1);
    yA=xA; xB=xA; yB=xA; xC=xA; yC=xA; vxA=xA; vyA=xA; vxB=xA; vyB=xA; vxC=xA;
    vyC=xA; axA=xA; ayA=xA; axB=xA; ayB=xA; axC=xA; ayC=xA;
    dth3=xA; dth4=xA; ddth3=xA; ddth4=xA;
    
    [K1, K2, K3, K4, K5] = getKBarRelations(r1, r2, r3, r4);
    
    for i=1:length(th2simul)
        th2 = th2simul(i);
        A = cos(th2)-K1-K2*cos(th2)+K3;
        B = -2*sin(th2);
        C = K1-(K2+1)*cos(th2)+K3;
        th4 = 2*atan((-B-sqrt(B^2-4*A*C))/(2*A));
        th3 = getTheta3(K1, K4, K5, th2);
        %----------------------------------
        % Local reference frame:
        xAr = r2*cos(th2);
        yAr = r2*sin(th2);
        xBr = xAr+r3*cos(th3);
        yBr = yAr+r3*sin(th3);
        xCr = r2*cos(th2)+rcx*cos(th3)-rcy*sin(th3);
        yCr = r2*sin(th2)+rcx*sin(th3)+rcy*cos(th3);
        %----------------------------------
        % Global reference frame:        
        Mtf = [cos(th0) -sin(th0); sin(th0) cos(th0)];
    
        AXY = Mtf*[xAr; yAr]+[x0;y0];
        BXY = Mtf*[xBr; yBr]+[x0;y0];
        CXY = Mtf*[xCr; yCr]+[x0;y0];
        
        xA(i) = AXY(1); yA(i)=AXY(2);
        xB(i) = BXY(1); yB(i)=BXY(2);
        xC(i) = CXY(1); yC(i)=CXY(2);
        
        if calcs==1
            syms omega3 omega4 alpha3 alpha4
            %%%
            % angular vel.
            w2 = [0 0 omega2];
            w3 = [0 0 omega3];
            w4 = [0 0 omega4];
            % angular acc.
            a2 = [0 0 alpha2];
            a3 = [0 0 alpha3];
            a4 = [0 0 alpha4];
    
            rO2A = [r2*cos(th2) r2*sin(th2) 0];
            rAB = [r3*cos(th3) r3*sin(th3) 0];
            rO4B = [r4*cos(th4) r4*sin(th4) 0];
            rAC = [rcx rcy 0];
            % Relative velocities equations
            vB = cross(w4, rO4B);
            vA = cross(w2, rO2A);
            vB_A = cross(w3, rAB);
            % Solve the system of equations
            sol1 = solve(vA+vB_A-vB==0,omega3, omega4);
    
            omega3 = eval(sol1.omega3);
            omega4 = eval(sol1.omega4);
            vB = eval(vB);
            vC = vA+cross(w3, rAC);
            vC = eval(vC);
    
            % Accelerations:
            aBn = cross(w4,vB);
            aBt = cross(a4,rO4B);
            aAn = cross(w2,vA);
            aAt = cross(a2,rO2A);
            aB_An = cross(w3,vB_A);
            aB_At = cross(a3,rAB);
    
            sol2 = solve(aAn+aAt+aB_An+aB_At-aBn-aBt==0,alpha3, alpha4);
        
            alpha3 = eval(sol2.alpha3);
            alpha4 = eval(sol2.alpha4);
            aC_An = cross(w3,cross(w3, rAC));
            aC_At = cross(a3,rAC);
            aCn = aAn+aC_An;
            aCt = aAt+aC_At;
            aCn = eval(aCn);
            aCt = eval(aCt);
            aC = aCn+aCt;

            aA = aAn+aAt;    
            aB = eval(aBn)+eval(aBt);    
            %----------------------------------
            % Global reference frame:
            vAXY = Mtf*vA(1:2)';
            vBXY = Mtf*vB(1:2)';
            vCXY = Mtf*vC(1:2)';
    
            aAXY = Mtf*aA(1:2)';
            aBXY = Mtf*aB(1:2)';
            aCXY = Mtf*aC(1:2)';
    
            vxA(i) = vAXY(1); vyA(i) = vAXY(2);
            vxB(i) = vBXY(1); vyB(i) = vBXY(2);
            vxC(i)=vCXY(1); vyC(i) = vCXY(2);
    
            axA(i) = aAXY(1); ayA(i) = aAXY(2);
            axB(i) = aBXY(1); ayB(i) = aBXY(2);
            axC(i) = aCXY(1); ayC(i) = aCXY(2);
    
            dth3(i) = omega3; dth4(i) = omega4;
            ddth3(i) = alpha3; ddth4(i) = alpha4;
        end
    end
    if calcs==1    
        TPs = table(xA,yA,xB,yB,xC,yC,vxA,vyA,vxB,vyB,vxC,vyC,axA,ayA,axB,ayB,axC,ayC,dth3,dth4,ddth3,ddth4);        
    else
        TPs = table(xA,yA,xB,yB,xC,yC);
    end
end
