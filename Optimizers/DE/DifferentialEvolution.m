function [bestPoint, bestFitness] = DifferentialEvolution(func, bounds, np, gen_max, iF, CR, useBest, y, fixCriteria)
    %DifferentialEvolution Implementation of Differential Evolution: https://www1.icsi.berkeley.edu/~storn/code.html
    %   This is DE/{rand_best}/1/bin with {none_generation_vector} Dither
    %   Inputs: func is a function pointer to the objective function
    %           bounds is a matrix of Dx2 size. Every i row is [L_i, U_i] => lower and upper bound of dimension i
    %           np is the number of individuals in the population
    %           gen_max is the number of generations to execute
    %           iF refers to the F factor. Set iF in (0, 2] for no dither, i.e., iF = F
    %                                      Set iF = -1 for per-generation dither, i.e., F_G in [0.5, 1.0]
    %                                      Set iF = -2 for per-vector dither, i.e., F_v in [0.5, 1.0]
    %           CR is the per-bit crossover-rate and must be in [0, 1]
    %           Set useBest to true in order to use DE/best and to false in order to use DE/rand
    %           y is the number of difference vectors used, usually 1 or 2.
    %           fixCriteria is a function handle to improve 10% of random solutions (leave it as [] to ignore)
    %   Output: The best point found and its value
    %   Main References: R. Storn & K. Proce. Differential evolution: a simple and efficient heuristic for global optimization over continuous spaces. Journal of global optimization, 11(4), 341-359, 1997.
    %                    https://www1.icsi.berkeley.edu/~storn/code.html ; 
    %                    https://pablormier.github.io/2017/09/05/a-tutorial-on-differential-evolution-with-python/
    %   Coded by N.C. Cruz (University of Almeria, Spain), 2020
    if iF>-1
       F = iF;
    end
    numR = 1+2*y;%Number of random individuals involved in each mutation (every vector results from two... +1 for the /rand flavour)
    
    nVars = size(bounds, 1);
    [pop, popValues, bestIndex] = initializePopulation(np, func, bounds, nVars, fixCriteria); % Every column is an individual. As many rows as variables / dimensions
    
    newPop = zeros(nVars, np); % x2 in the paper
    newPopValues = zeros(1, np);
   
    for count=1:1:gen_max % Number of cycles
        if iF== -1 % per-generation dither
            F = rand()*0.5 + 0.5;  % F from the interval [0.5, 1.0] randomly for each generation
        end
        newBestIndex = -1;
        newBestFitness = inf;
        for i=1:1:np % Iterate through the whole population
            tickets = 1:1:np;
            tickets(i) = []; % Any chosen vector must be different from each other and the running index
            r = tickets(randperm(np-1, numR)); % The every difference vector needs 2 defining ones, and 1+ comes from the base one (for rand mode)
            if useBest
                r(1) = bestIndex;
            end
            k = randi(nVars); % Randomly pick the parameter to be inherited always
            
            trial = pop(:, i);
            if iF== -2 % per-vector dither
                F = rand()*0.5 + 0.5;  % F from the interval [0.5, 1.0] randomly for each generation
            end
            for j=1:1:nVars
                if rand()<=CR || j==k
                    trial(j) = pop(j,r(1)) + F*( sum(pop(j, r(2:2:numR)) - pop(j, r(3:2:numR))) ); % Mutant vector 'just in time'; 1V: = pop(j,r(1)) + F*(pop(j,r(2)) - pop(j,r(3)));
                    if trial(j) < bounds(j,1) % Fix limits!!
                        trial(j) = bounds(j,1);
                    elseif trial(j) > bounds(j,2)
                        trial(j) = bounds(j,2);
                    end
                end % There is no need to use 'else' as trial is equal to pop by default
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
    pop = zeros(nVars, np);
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
