function J = light_obj_func(x, context)
    XY = get_XY(x, context);
    if context.prescribed == 3
        h2 = Sequence(x(7:6+context.lenRefXY));
    else
        h2 = 0;
    end
    
    if context.costType == 1 % NSDV:
        A = GetShapeVector(XY(:,1),XY(:,2),1);
        % Scale checking!!
        L_opt = normalized_L(XY);
        X_sca = reshape(x*(context.L_des/L_opt), length(x), 1); % Ensuring that you get a column (TLBO for instance stores x as a row)
        if any(X_sca(1:6) < context.barBounds(1:6, 1)) || any(X_sca(1:6) > context.barBounds(1:6, 2))
            discard = 10000*1;
        else
            discard = 0;
        end
        J=norm(A - context.refShapeVector) + h2*1000000 + discard;
    else % Absolute:
        Ji=zeros(context.lenRefXY, 1);
        for i=1:context.lenRefXY
            Ji(i)=(XY(i,1) - context.XY_des(i, 1))^2 + (XY(i,2) - context.XY_des(i, 2))^2;
        end
        J=sum(Ji) + h2*1000000;
    end
end
