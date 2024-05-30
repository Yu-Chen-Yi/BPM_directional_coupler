function [bestSolution, bestFitness] = greyWolfOptimization(fitnessFunction, numVariables, lowerBound, upperBound, maxGenerations, populationSize)
    % initial groups
    alphaWolf = struct('position', [], 'fitness', inf);
    betaWolf = struct('position', [], 'fitness', inf);
    deltaWolf = struct('position', [], 'fitness', inf);
    wolves = repmat(struct('position', [], 'fitness', inf), 1, populationSize);
    
    % initial wolfs
    for i = 1:populationSize
        wolves(i).position = rand(1, numVariables) .* (upperBound - lowerBound) + lowerBound;
        wolves(i).fitness = fitnessFunction(wolves(i).position);
        
        % update alpha, beta and delta wolfs
        if wolves(i).fitness < alphaWolf.fitness
            deltaWolf = betaWolf;
            betaWolf = alphaWolf;
            alphaWolf = wolves(i);
        elseif wolves(i).fitness < betaWolf.fitness
            deltaWolf = betaWolf;
            betaWolf = wolves(i);
        elseif wolves(i).fitness < deltaWolf.fitness
            deltaWolf = wolves(i);
        end
    end
    
    for generation = 1:maxGenerations
        a = 2 - generation * (2 / maxGenerations);
        
        for i = 1:populationSize
            for j = 1:numVariables
                A1 = 2 * a * rand() - a;
                C1 = 2 * rand();
                D_alpha = abs(C1 * alphaWolf.position(j) - wolves(i).position(j));
                X1 = alphaWolf.position(j) - A1 * D_alpha;
                
                A2 = 2 * a * rand() - a;
                C2 = 2 * rand();
                D_beta = abs(C2 * betaWolf.position(j) - wolves(i).position(j));
                X2 = betaWolf.position(j) - A2 * D_beta;
                
                A3 = 2 * a * rand() - a;
                C3 = 2 * rand();
                D_delta = abs(C3 * deltaWolf.position(j) - wolves(i).position(j));
                X3 = deltaWolf.position(j) - A3 * D_delta;
                
                wolves(i).position(j) = (X1 + X2 + X3) / 3;
                
                % limit the wolfs in the boundary condition
                wolves(i).position(j) = max(wolves(i).position(j), lowerBound(j));
                wolves(i).position(j) = min(wolves(i).position(j), upperBound(j));
            end
            
            % uadata wolfs' fitness
            wolves(i).fitness = fitnessFunction(wolves(i).position);
            
            % update alpha, beta and delta wolfs
            if wolves(i).fitness < alphaWolf.fitness
                deltaWolf = betaWolf;
                betaWolf = alphaWolf;
                alphaWolf = wolves(i);
            elseif wolves(i).fitness < betaWolf.fitness
                deltaWolf = betaWolf;
                betaWolf = wolves(i);
            elseif wolves(i).fitness < deltaWolf.fitness
                deltaWolf = wolves(i);
            end
        end
        
    end
    
    % retrun the bestsolution and bestfitness
    bestSolution = alphaWolf.position;
    bestFitness = alphaWolf.fitness;
end