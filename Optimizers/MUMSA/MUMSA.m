function [bestPoint, bestFitness] = MUMSA(func, bounds, np, itermax, f, cp, mp, range, fixCriteria)
    %MUMSA Implementation of the MUMSA (Malaga University Mechanism Synthesis Algorithm) method
    %   The MUMSA (Malaga Univ. Mechanism Synthesis Algorithm) was proposed in
    %   J.A. Cabrera, A. Ortiz, F. Nadal & J.J. Castillo. An evolutionary
    %   algorithm for path synthesis of mechanisms. Mechanism and Machine
    %   Theory, 46(2), 127-141, 2011. Up to me, it is mainly Differential
    %   Evolution in its flavour DE/best/1/bin without dither and plus a new mutation component
    %   Inputs: func is a function pointer to the objective function
    %           bounds is a matrix of Dx2 size. Every i row is [L_i, U_i] => lower and upper bound of dimension i
    %           np is the number of individuals in the population
    %           itermax is the number of generations to execute
    %           f, in R range (0, 2], controls the amplification of the differential variation
    %           cp is the per-bit crossover-rate and must be in [0, 1]
    %           mp defines the mutation probability of adding disturbances
    %           in +/-range (which must be entered in absolute value)
    %           fixCriteria is a function handle to improve 10% of random solutions (leave it as [] to ignore)
    %   Output: The best point found and its value
    %   Coded by N.C. Cruz (University of Almeria, Spain), 2020
    nVars = size(bounds, 1); % The number of rows should be equal to the number of variables
    
    [pop, popValues, bestIndex] = initializePopulation(np, func, bounds, nVars, fixCriteria); % Every column is an individual. As many rows as variables / dimensions
    
    newPop = zeros(nVars, np); % According to Differential Evolution, new individuals do not affect until the next generation
    newPopValues = zeros(1, np);
    
    for count=1:1:itermax
        newBestIndex = -1; % To keep track of the new best individual without looking for it
        newBestFitness = inf;
        for i=1:1:np
            tickets = 1:1:np;
            tickets(i) = []; % r1 and r2 are chosen to be different from each other and the running index
            r = tickets(randperm(np-1, 2));
            k = randi(nVars); % Paramter forced to be inherited 
            
            trial = pop(:, i);
            for j=1:1:nVars
                if rand()<=cp || j==k
                    trial(j) = pop(j, bestIndex) + f*(pop(j,r(1)) - pop(j,r(2))); % Using the mutant vector 'just in time'
                end
                if rand()<mp % MUMSA Mutation of the new component (either from the mutant vector or from the population)
                    my_min = trial(j)-range;
                    my_max = trial(j)+range;
                    trial(j) = rand()*(my_max-my_min)+my_min; % choosing a real value in range (x_i, +-range)
                end
                if trial(j) < bounds(j,1) % Fix limits in any case!!
                    trial(j) = bounds(j,1);
                elseif trial(j) > bounds(j,2)
                    trial(j) = bounds(j,2);
                end
            end
            
            score = func(trial);
            if score<=popValues(i) % It improves the initial value
                newPop(:,i) = trial;
                newPopValues(i) = score;
            else                    % It does not... keep the original
                newPop(:,i) = pop(:, i);
                newPopValues(i) = popValues(i);
            end
            if newPopValues(i)<newBestFitness % Keeping the best of the generation up to date anyway
                newBestIndex = i;
                newBestFitness = newPopValues(i);
            end
        end
        bestIndex = newBestIndex;
        [newPop, pop] = swap(newPop, pop);
        [newPopValues, popValues] = swap(newPopValues, popValues);
    end
    
    bestPoint = pop(:,bestIndex);
    bestFitness = popValues(bestIndex);
end

% Internal auxiliary functions:

function [pop, popValues, bestIndex] = initializePopulation(np, func, bounds, nVars, fixCriteria)
    pop = zeros(nVars, np); % Every column is an individual with as many rows as variables
    popValues = zeros(1, np);
    
    bestIndex = 1;
    bestVal = inf;
    
    for i=1:1:np
        new_point = rand(nVars,1).*(bounds(:,2)-bounds(:,1)) + bounds(:,1); % [0,1]*(max-min) + min
        if ~isempty(fixCriteria) && mod(i, 10)==0
            new_point = fixCriteria(new_point);
        end
        pop(:, i) = new_point;
        popValues(i) = func(new_point);
        if popValues(i)<bestVal
            bestIndex = i;
            bestVal = popValues(i);
        end
    end
end

function [b, a] = swap(a, b)
    % See: https://es.mathworks.com/matlabcentral/answers/362680-how-to-swap-values-of-two-variables
end
