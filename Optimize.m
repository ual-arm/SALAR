function [context, resultPack] = Optimize(desiredPath, th2_n, bounds, costType, probType, optimConfig)
    resultPack = [];
    context = BuildOptimContext(desiredPath, th2_n, costType, probType, bounds);
    
    optimNames = LoadOptimNames();
    [bounds, optimConfig] = Sanitize(bounds, optimConfig, costType, probType, th2_n, optimNames);
    
    if isfield(optimConfig, optimNames.OptimizationFMCName)
        fmcConfig = optimConfig.(optimNames.OptimizationFMCName);
        resultPack.(optimNames.fmcResultsName) = launchFMinCon(context, bounds, fmcConfig, optimNames);
    end
    if isfield(optimConfig, optimNames.OptimizationTLBOName)
        tlboConfig = optimConfig.(optimNames.OptimizationTLBOName);
        resultPack.(optimNames.tlboResultsName) = launchExternalOptimizer(context, bounds, 'TLBO', ...
            'TLBOSearch', tlboConfig.(optimNames.OptimizationTLBORestartsName), tlboConfig, optimNames);
    end
    if isfield(optimConfig, optimNames.OptimizationMUMSAName)
        mumsaConfig = optimConfig.(optimNames.OptimizationMUMSAName);
        resultPack.(optimNames.mumsaResultsName) = launchExternalOptimizer(context, bounds, 'MUMSA', ...
            'MUMSA', mumsaConfig.(optimNames.OptimizationMUMSARestartsName), mumsaConfig, optimNames);
    end
    if isfield(optimConfig, optimNames.OptimizationDEName)
        deConfig = optimConfig.(optimNames.OptimizationDEName);
        resultPack.(optimNames.deResultsName) = launchExternalOptimizer(context, bounds, 'DE', ...
            'DE', deConfig.(optimNames.OptimizationDERestartsName), deConfig, optimNames);
    end
end

% Internal auxiliary functions:
function context = BuildOptimContext(desiredPath, th2_n, costType, probType, bounds)
    context.XY_des = desiredPath;
    context.lenRefXY = size(desiredPath, 1); % The number of rows is equal to the number of points
    context.prescribed = probType;
    context.th2_n = th2_n;
    context.costType = costType;
    
    context.L_des = normalized_L([desiredPath(:, 1), desiredPath(:, 2)]);
    context.isRefClockWise = IsClockWiselyOrdered(reshape(context.XY_des', 1, 2*size(context.XY_des, 1)));
    context.barBounds = bounds(1:6,:);
    context.refShapeVector = GetShapeVector(desiredPath(:,1), desiredPath(:,2), true);
    if costType == 1 % NSDV-Way
        context.errorThreshold = sqrt(2*context.lenRefXY - 1);
    else
        context.errorThreshold = 2*context.lenRefXY*(120^2);
    end
end

function [bounds, optimConfig] = Sanitize(bounds, optimConfig, costType, probType, th2_n, names)
    if costType==1
        if probType==1
            bounds=bounds(1:6,:);
            if IncludeX0(optimConfig)
                optimConfig.(names.OptimizationFMCName).(names.OptimizationFMCX0Name) = ...
                    optimConfig.(names.OptimizationFMCName).(names.OptimizationFMCX0Name)(1:6);
            end
        elseif probType==2
            bounds=bounds(1:7,:);
            if IncludeX0(optimConfig)
                optimConfig.(names.OptimizationFMCName).(names.OptimizationFMCX0Name) = ...
                    optimConfig.(names.OptimizationFMCName).(names.OptimizationFMCX0Name)(1:7);
            end
        elseif probType==3
            bounds=bounds(1:6+length(th2_n),:);
            if IncludeX0(optimConfig)
                optimConfig.(names.OptimizationFMCName).(names.OptimizationFMCX0Name) = ...
                    optimConfig.(names.OptimizationFMCName).(names.OptimizationFMCX0Name)(1:6+length(th2_n));
            end
        end
    end
    % Nested function that checks if the optimConfig has an initial point for FMC
    function hasX0 = IncludeX0(optimConfig)
        hasX0 = false;
        if isfield(optimConfig, names.OptimizationFMCName) && isfield(optimConfig.(names.OptimizationFMCName), names.OptimizationFMCX0Name)
            hasX0 = true;
        end
    end
end

function solution = FixSolution(solution, prescribed_mode, np)
    % Make the second bar be shortest:
    idShortest = find(solution(1:4) == min(solution(1:4)), 1);
    if idShortest~=2
        swap = solution(2);
        solution(2) = solution(idShortest);
        solution(idShortest) = swap;
    end
    % Let's order the angles, if needed:
    if prescribed_mode == 3
        th2i = solution(7:6+np);
        sorted_angles = sort(th2i);
        solution(7:6+np) = sorted_angles;
    end
end

function emptyResult = BuildOptimResult(nameSet)
    emptyResult.(nameSet.bestRawResultName) = [];
    emptyResult.(nameSet.bestRawFuncValName) = inf;
    emptyResult.(nameSet.meanQualityName) = 0;
    emptyResult.(nameSet.stdQualityName) = 0;
    emptyResult.(nameSet.meanTimeName) = 0;
end

function fmcResults = launchFMinCon(context, bounds, fmcConfig, nameSet)
    fmcResults = BuildOptimResult(nameSet);
    
    tamVecs = size(bounds(:,1));
    vecQuality = zeros(1, fmcConfig.(nameSet.OptimizationFMCRestartsName));
    vecTime = zeros(1, fmcConfig.(nameSet.OptimizationFMCRestartsName));
    Lb = bounds(:, 1);
    Ub = bounds(:, 2);
    
    for i=1:1:fmcConfig.(nameSet.OptimizationFMCRestartsName)
        warning('off'); % To avoid internal messages from FMinCon (Specific warning ID seems buggy)
        tA = tic;
        if i==1 && fmcConfig.(nameSet.OptimizationFMCUStartName)
            x0 = fmcConfig.(nameSet.OptimizationFMCX0Name);
            A = zeros(3,length(x0));            
            A(:,1:4) = [-1 1 -1 1; -1 1 1 -1; 1 1 -1 -1];            
            b = ones(length(A(:,1)),1)*0.00001; % To achieve < instead of <=
            options = optimoptions(@fmincon,'Display', 'iter', 'Algorithm', 'interior-point', 'MaxFunctionEvaluations', 1e7);
            [x, fval] = fmincon(@(x)light_obj_func(x, context), x0, A, b, [], [], Lb, Ub, [], options);
        else
            x0 = FixSolution(rand(tamVecs).*(Ub - Lb) + Lb, context.prescribed, context.lenRefXY);
            [x, fval] = fmincon(@(x)armored_obj_func(x, context), x0, [], [], [], [], Lb, Ub);
        end
        tB = toc(tA);
        warning('on');
        vecQuality(i) = fval;
        vecTime(i) = tB;
        if fval < fmcResults.(nameSet.bestRawFuncValName)
            fmcResults.(nameSet.bestRawResultName) = x;
            fmcResults.(nameSet.bestRawFuncValName) = fval;
        end
    end
    
    fmcResults.(nameSet.meanQualityName) = mean(vecQuality);
    fmcResults.(nameSet.stdQualityName) = std(vecQuality);
    fmcResults.(nameSet.meanTimeName) = mean(vecTime);
end

function call = BuildExternOptimCall(context, bounds, methodName, config, nameSet)
    switch methodName
        case 'TLBO'
            call = @() TLBOSearch(@(x)armored_obj_func(x, context), config.(nameSet.OptimizationTLBOPopSizeName), ...
            config.(nameSet.OptimizationTLBONCyclesName), bounds, @(x)FixSolution(x, context.prescribed, context.lenRefXY)); 
        case 'MUMSA'
            call = @() MUMSA(@(x)armored_obj_func(x, context), bounds, config.(nameSet.OptimizationMUMSANPName), ...
            config.(nameSet.OptimizationMUMSAItermaxName), config.(nameSet.OptimizationMUMSAFName), config.(nameSet.OptimizationMUMSACPName),...
            config.(nameSet.OptimizationMUMSAMPName), config.(nameSet.OptimizationMUMSARangeName), @(x)FixSolution(x, context.prescribed, context.lenRefXY));
        case 'DE'
            call = @() DifferentialEvolution(@(x)armored_obj_func(x, context), bounds, config.(nameSet.OptimizationDENPName), ...
            config.(nameSet.OptimizationDEGenMaxName), config.(nameSet.OptimizationDEFName), config.(nameSet.OptimizationDECRName), ...
            config.(nameSet.OptimizationDEXName), config.(nameSet.OptimizationDEYName), @(x)FixSolution(x, context.prescribed, context.lenRefXY));
    end
end

function results = launchExternalOptimizer(context, bounds, methodName, folderName, numRestarts, config, nameSet)
    newPath = ['.', filesep, 'Optimizers', filesep, folderName];
    addpath(newPath);
    
    results = BuildOptimResult(nameSet);
    vecQuality = zeros(1, numRestarts);
    vecTime = zeros(1, numRestarts);
    
    call = BuildExternOptimCall(context, bounds, methodName, config, nameSet);
    
    for i=1:1:numRestarts
        tA = tic;
        [x, fval] = call();
        vecTime(i) = toc(tA);
        vecQuality(i) = fval;
        if fval < results.(nameSet.bestRawFuncValName)
            results.(nameSet.bestRawResultName) = x;
            results.(nameSet.bestRawFuncValName) = fval;
        end
    end
    
    if size(results.(nameSet.bestRawResultName), 2) > 1 % Only column vectors are admitted as a result, but some methods may return rows
        results.(nameSet.bestRawResultName) = results.(nameSet.bestRawResultName)';
    end
    
    results.(nameSet.meanQualityName) = mean(vecQuality);
    results.(nameSet.stdQualityName) = std(vecQuality);
    results.(nameSet.meanTimeName) = mean(vecTime);
    
    rmpath(newPath);
end
