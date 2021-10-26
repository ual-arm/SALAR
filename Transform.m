function [X1_mod, th2_prescrib_mod] = Transform(X1, context, isClosed)
    np = context.lenRefXY;
    XY_des = context.XY_des;
    error = inf;
    XY1 = get_XY(X1, context);
    switch context.prescribed
        case 1
            th2_n=[];
        case 2
            th2_n=X1(7);
        case 3
            th2_n=X1(7:6+np);
    end

    L_opt=normalized_L(XY1);
    %================ Step 1: Scale
    x2=[ X1(1:6)*(context.L_des/L_opt); th2_n; 0; 0; 0 ];
    % Centroid of the desired curve:
    OXdes = sum(XY_des(:,1))/length(XY_des);
    OYdes = sum(XY_des(:,2))/length(XY_des);
    th0_i = 0;
    angulos_original_t12 = context.th2_n;
    angulos_original_t3 = th2_n;
    
    th2_prescrib_mod = context.th2_n;
    for i=1:360
        %================ Step 2: Rotate (x3 candidate)
        x3 = [ x2(1:6); th2_n; 0; 0; th0_i ];
        XY3 = get_XY(x3, context);
        OXgen = sum(XY3(:,1))/length(XY3);
        OYgen = sum(XY3(:,2))/length(XY3);
        %===================== Step 3: Translate (x4 candidate -> evaluate each th2_1
        x4 = [x3(1:6); th2_n; (OXdes-OXgen); (OYdes-OYgen); th0_i ];
        if isClosed && context.prescribed == 3
            for j=1:np
                % We try every possible starting point (th2_1) and choose that with less Absolute Error, for each th0:
                x4(7:6+np) = circshift(angulos_original_t3,j);
                % Calculate the absolute error:
                XY4 = get_XY(x4, context);
                error_ij = sum( ( XY4(:,1) - XY_des(:,1) ).^2 + ( XY4(:,2) - XY_des(:,2) ).^2 );
                if error_ij < error
                    error = error_ij;
                    th0_out = th0_i;
                    th2_n_out = x4(7:6+np); % This is the one which minimizes AE at the optimal th0
                end
            end
        elseif isClosed
            th2_n_out = th2_n; % Close but not of type 3 => Acting on prescribed angles only
            for j=1:np
                % We try every possible starting point (th2_1) and choose that with less Absolute Error, for each th0:
                Buffer_th2_prescrib = circshift(angulos_original_t12,j);
                context.th2_n = Buffer_th2_prescrib; % Make effective the change
                XY4 = get_XY(x4, context);
                error_ij = sum( ( XY4(:,1) - XY_des(:,1) ).^2 + ( XY4(:,2) - XY_des(:,2) ).^2 );
                if error_ij < error
                    error = error_ij;
                    th2_prescrib_mod = Buffer_th2_prescrib;
                    th0_out = th0_i;
                end
            end
        else % Neither closed nor of type 3: Do not modify th_n or the prescribed angles
            th2_n_out = th2_n;            
            XY4 = get_XY(x4,context);
            error_ij = sum( ( XY4(:,1) - XY_des(:,1) ).^2 + ( XY4(:,2) - XY_des(:,2) ).^2 );
            if error_ij < error
                th0_out = th0_i;
                error = error_ij;
            end
        end
        th0_i = th0_i + deg2rad(1);
    end
    %=========== Step 4: Scaling + rotation + translation for the optimal th21.
    if isClosed && context.prescribed == 2
        context.th2_n = th2_prescrib_mod;
    end
    x3 =[ x2(1:6); th2_n_out; 0; 0; th0_out ];
    XY3 = get_XY(x3, context);
    OXgen = sum(XY3(:,1))/length(XY3);
    OYgen = sum(XY3(:,2))/length(XY3);
    X1_mod=[x3(1:6); th2_n_out; (OXdes-OXgen); (OYdes-OYgen); th0_out];
end
