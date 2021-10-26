function J = armored_obj_func(x, context)
    penaltyFactor = context.errorThreshold;
    
    bars = x(1:4);
    
    id_short = find(bars==min(bars)); % Look for every instance! Not ,1...
    if ismember(2, id_short)
        id_short = 2; % There might be several with the shortest distance, but we don't care as long as R2 is one of them
    else
        id_short = id_short(1);
    end
    short = bars(id_short);
    desiredShort = bars(2); % We will compare the previous val with the desired.
    
    id_long = find(bars==max(bars), 1);
    long = bars(id_long);
    
    bars([id_short, id_long]) = [];
    sum_meds = sum(bars);
    
    if std(x(1:4))==0 || ( id_short==2 && (short+long)<=sum_meds ) % If all the bars are the same OR R2 is the shortest and you meet Grashof
        XY=get_XY(x, context);
        if sum(isnan(XY(:)))==0 % If the simulator was able to run flawlessly... Ok, go on
            seqMatters = 0;
            if context.prescribed==3
                seqMatters = Sequence(x(7:6+context.lenRefXY));
            end
            if context.costType == 1
                L_opt = normalized_L(XY);
                X_sca = reshape(x*(context.L_des/L_opt), length(x), 1); % Ensuring you get a column vector to compare (TLBO uses rows for instance)
                badScale = 0;
                if any(X_sca(1:6) < context.barBounds(1:6, 1)) || any(X_sca(1:6) > context.barBounds(1:6, 2))
                    badScale = penaltyFactor;
                end
                badSense = 0;
                isMineClockWise = IsClockWiselyOrdered(reshape(XY', 1, 2*size(XY, 1)));
                if context.isRefClockWise ~= isMineClockWise
                    badSense = penaltyFactor;
                end
                A = GetShapeVector(XY(:,1),XY(:,2),1);
                J = norm(A-context.refShapeVector) + penaltyFactor*(seqMatters + badScale + badSense);
                if isnan(J)
                    J = penaltyFactor^4; % VERY Misteriously bad solution
                end
            elseif context.costType == 0
                J = sum( (XY(:,1)-context.XY_des(:,1)).^2 + (XY(:,2)-context.XY_des(:,2)).^2 ) + penaltyFactor*seqMatters;
            else
                error('Unknown cost method');
            end
        else
            J = penaltyFactor^3; % Misteriously bad solution
        end
    else
        J = penaltyFactor;
        if desiredShort > short
           J = J + penaltyFactor*(desiredShort-short)/desiredShort; 
        end
        if std(x(1:4))~=0 && id_long == 2 % If not all are the same and the longest should not be there
            J = J + penaltyFactor;
        end
        if std(x(1:4))~=0 && (desiredShort+long) > sum_meds && id_short == 2 % The last term is to ensure that sum_meds is appropriate
           J = J + penaltyFactor*( (desiredShort+long) - sum_meds )/(desiredShort+long); 
        end
    end
end
