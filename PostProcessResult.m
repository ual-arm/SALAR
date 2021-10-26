function resultPack = PostProcessResult(optimContext, optimResult, isClosed)
    nameSet = LoadOptimNames();
    
    [X1, bestRawVal, log] = ScanOptimResult(optimResult, nameSet);
    resultPack.(nameSet.bestXRawValName) = bestRawVal;
    resultPack.(nameSet.LogName) = log; % Storing the log
    
    if bestRawVal < optimContext.errorThreshold % If the result deserves to be studied:
        th2_prescrib = optimContext.th2_n; % (If absoute, current_th2 does not vary)
        
        if optimContext.costType == 1 % NSDV 
            [X1, th2_prescrib] = Transform(X1, optimContext, isClosed);
            optimContext.th2_n = th2_prescrib; % Update the context for local simulation if needed
        end
    
        if optimContext.prescribed == 3
            th2_prescrib = X1(7:6+optimContext.lenRefXY);
        end
    
        XYachieved = get_XY(X1, optimContext);

        resultPack.(nameSet.bestXName) = X1;
        resultPack.(nameSet.XYachievedName) = XYachieved;
        resultPack.(nameSet.UpdatedTh2Name) = th2_prescrib;
    
        resultPack.(nameSet.costTDName) = ComputeTDError(XYachieved, optimContext.XY_des);
        optimContext.costType = 0; % AE Cost
        resultPack.(nameSet.costAEName) = armored_obj_func(X1, optimContext);
        optimContext.costType = 1; % NSD-V cost
        resultPack.(nameSet.costNSDVName) = armored_obj_func(X1, optimContext);
    end
end

% Internal auxiliary functions:
function [X1, bestVal, log] = ScanOptimResult(optimResult, nameSet)
    log = '';
    X1 = [];
    winPos = 0;
    bestVal = inf;
    
    optNames = {'FMinCon', 'TLBO', 'MUMSA', 'DE'};
    resultNames = {nameSet.fmcResultsName, nameSet.tlboResultsName, nameSet.mumsaResultsName, nameSet.deResultsName}; 
    formatString = '\n----------------\n%s:\nBest Solution: %s\nValue of the best solution: %f\nMean value: %f (STD: %f)\nMean time: %f (s)\n';
    
    for i=1:1:length(resultNames)
        if isfield(optimResult, resultNames{i})
            strc = optimResult.(resultNames{i});
            fragment = sprintf(formatString, optNames{i}, mat2str(strc.(nameSet.bestRawResultName)'), strc.(nameSet.bestRawFuncValName), ... 
                strc.(nameSet.meanQualityName), strc.(nameSet.stdQualityName), strc.(nameSet.meanTimeName));
            log = strcat(log, fragment);
            
            bestValCandidate = strc.(nameSet.bestRawFuncValName);
            if bestValCandidate < bestVal
                bestVal = bestValCandidate;
                X1 = strc.(nameSet.bestRawResultName);
                winPos = i;
            end
        end
    end
    
    if winPos>0
        log = strcat(log, sprintf('\n\nWINNER: %s\n', optNames{winPos}));
    else
        log = strcat(log, sprintf('\n\nWINNER: <NONE>\n'));
    end
end

function cost_TD = ComputeTDError(XYachieved, XYdesired) % The Nadal et al's cost, but with the Matlab's built-in Turning Functions
    try
        AA=polyshape(XYachieved(:, 1), XYachieved(:, 2));
        BB=polyshape(XYdesired(:, 1), XYdesired(:, 2));
        cost_TD=turningdist(AA, BB);
    catch
        cost_TD=Inf;
    end
end
