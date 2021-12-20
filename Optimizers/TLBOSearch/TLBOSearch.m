% Teaching-Learning-Based Optimization algorithm compact implementation for unconstrained optimization problems for Matlab
% This algorithm was proposed by Rao et al (2012).
% INPUTS:   # A pointer to the objective function to be used (pay attention 
%           to any previuos/later initialization of it, e.g. specific MEX functions...)
%           # The size of the population, the number of generations to run,
%           the number of optimization variables M an a Mx2 matrix with the
%           limits of the searchspace (lower and upper bound at every row)
% OUTPUTS:  # The vector of the best solution found along the process and
%           its value according to the given objective function
% OPEN-SOURCE, FREE TO USE, SHARE AND MODIFY. Implemented by N.C. Cruz, University of Almeria, Spain.

function [bestSolution, FEVAL] = TLBOSearch(readyCostFunction, populationSize, numGenerations, bounds, fixCriteria)% _STEP 1_ Algorithm configuration
    
    numVariables = size(bounds, 1);% The number of rows of the bounds should be equal to the number of variables
    
    % Initialize the population: _(STEP 2)_
    population = zeros(populationSize, numVariables); % Pre-allocating the room for populationSize individuals (rows) with numVariables subjects (cols)
    popCosts = zeros(populationSize, 1); % The value of the 'i' individual (a certain row at the population matrix)
    
    for i=1:1:populationSize % Loop for writing every indidividual (Row loop)
        for j=1:1:numVariables % Loop for writing every particular subject according to the limits (Column loop)
            population(i, j) = (bounds(j, 2)-bounds(j, 1))*rand(1,1) + bounds(j, 1);%(UBi-LBi)*(0,1) + LBi of that variable
        end
        if ~isempty(fixCriteria) && mod(i, 10)==0
            population(i, :) = fixCriteria(population(i, :));
        end
        
        popCosts(i) = readyCostFunction(population(i,:)); % Evaluating and saving the value of the objective function        
    end
    
    % Main loop _STEP 3_ STEP 4_ STE 5_ -> _STEP_3...
    for i=1:1:numGenerations
        % _STEP 3_: Teacher phase
        % Computing the population mean at every subject (column wise)
        M_D = mean(population); % For matrices, mean(X) is a row vector containing the mean value of each column.
        
        % Looking for the teacher (the best solution) as we do not keep the population ordered:
        [teacherFocus, ~] = LookForTheBest(populationSize, popCosts);
        
        M_new_D = population(teacherFocus, :);% The row is the individual
        r = rand(1, numVariables);% There is a random 'r_j' per subject (design variable)
        T_F = randi([1, 2], 1);% There is a global Teaching Factor at every generation
        Difference_D = r.*(M_new_D - T_F*M_D);
        
        for j =1:1:populationSize
            newSolution = population(j,:) + Difference_D;% Teaching...            
            newSolution = FixLimits(newSolution, bounds); % Keeping every design variable in its valid limits           
            newCost = readyCostFunction(newSolution);% Computing the value of the new solution
            
            if newCost < popCosts(j) % The new one is better, let's save it  
                population(j,:) = newSolution;
                popCosts(j) = newCost; %Update the vector of costs
            end % (It will be discarded otherwise)
        end
        % _STEP 4_: Learner/Student phase
        r = rand(1, numVariables);% There is a random 'r_j' per subject (design variable) <Common for all the stage>
        popSnapshot = population; % Copy: All interactions are before updating
        costSnapShot = popCosts; % Copy: All interactions are before updating
        for j = 1:1:populationSize
            validSelection = false;
            while ~validSelection % Seek until a different index is found (no self-interaction allowed)
                focus = randi([1, populationSize], 1);
                if focus ~= j
                    validSelection = true;
                end
            end
            if costSnapShot(j) < costSnapShot(focus)% The learner 'j' is better -> knowledge from 'j' to 'focus'
                newSolution = popSnapshot(j, :) + ( r .* (popSnapshot(j, :) - popSnapshot(focus, :)) );
            else % The learner 'focus' is better -> knowledge from 'focus' to 'j'
                newSolution = popSnapshot(j, :) + ( r .* (popSnapshot(focus, :) - popSnapshot(j, :)) );
            end
            newSolution = FixLimits(newSolution, bounds); % Keeping every design variable in its valid limits 
            newCost = readyCostFunction(newSolution);% Computing the value of the new solution
            
            if newCost < costSnapShot(j) % The new one is better, let's save it  
                population(j,:) = newSolution; % Altering the population matrix
                popCosts(j) = newCost; % Altering the vector of costs
            end % (It will be discarded otherwise)            
        end
        [population, popCosts] = RemoveDuplicates(readyCostFunction, population, popCosts, bounds); % Looking for a fully heterogeneous population
    end   
    % The procedure is over: lets pick-up the best individual:
    [bestFocus, FEVAL] = LookForTheBest(populationSize, popCosts);%It returns both the solution vector and its value according to the objective function
    bestSolution = population(bestFocus, :);
end

% -----------Internal Functions------:

function [bestFocus, bestCost] = LookForTheBest(populationSize, popCosts)
    bestFocus = 1;
    bestCost = popCosts(1);
    for j=2:1:populationSize
        if popCosts(j) < bestCost % This is a better teacher: let's update it
            bestFocus = j;
            bestCost = popCosts(j);
        end
    end
end

function solution = FixLimits(solution, bounds)
    for k=1:1:size(bounds, 1) % The number of rows of the vector should be the number of variables
        if solution(k) < bounds(k, 1) % Less than the minimum
            solution(k) =  bounds(k, 1);
        elseif solution(k) > bounds(k, 2) % More than the maximum
            solution(k) =  bounds(k, 2);
        end
    end
end

function [population, popCosts] = RemoveDuplicates(readyCostFunction, population, popCosts, bounds)%Mutating the individuals in case of duplicities
    numVariables = size(bounds, 1);
    
    for i = 1:1:size(population, 1)% The number of individuals is equal to the number of rows
        for j = i+1:1:size(population, 1)
            if isequal(population(i,:), population(j,:)) % The student at 'j' is equal to the student at 'i', let's change the 'j' one
                mutatedSolution = population(j, :);
                isNew = false;
                while ~isNew
                    subject = randi([1 numVariables], 1);%Selecting the subject (design variable) to be altered
                    mutatedSolution(subject) = rand()*( bounds(subject, 2) - bounds(subject, 1)) + bounds(subject, 1); % [0,1] * (max - min) + min
                    if ~ismember(mutatedSolution, population, 'rows')% Is it a real new individual?
                        isNew = true;
                    end
                end % Updating the population:
                population(j, :) = mutatedSolution;
                popCosts(j) = readyCostFunction(mutatedSolution);
            end
        end
    end
end
