function th3 = getTheta3(K1, K4, K5, th2)
    D = cos(th2) - K1 + K4*cos(th2) + K5;
    E = -2*sin(th2);
    F = K1 + (K4-1)*cos(th2) + K5;   
    th3 = 2*atan((-E-sqrt(E^2-4*D*F))/(2*D));
    % We take the open solution... (for cross solution take the +sqrt)
end
