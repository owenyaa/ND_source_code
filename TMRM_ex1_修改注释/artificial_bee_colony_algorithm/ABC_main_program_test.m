%% Author : GUANG_LIU  * owenyaa@gmail.com *
% Created Time : 2020-11-12 21:32
% Last Revised : GUANG_LIU ,2020-11-12
% Remark : The structure, function and variables of this program are explained:%
%          This program is to calculate the example 1 in the paper "A New Semi-%
%          analytical technique for Nonlinear Systems Based on Response Sensitivity%
%           Analysis", that is, the semi-analytical solution of equation %
%          \begin{equation}%
%          	\label{eq4.2}%
%          	\begin{cases}%
%          		\begin{aligned}%
%          		  & \dot{S}(t)+S^2(t)=1, t\geq 0 \\%
%          		  & S(0)=0                       %
%          		\end{aligned}%
%          	\end{cases}%
%          \end{equation}. The structure of the program is as follows:%
%          The basis function in this program is exponential function \chi_i(t)=e^{-it}%
%          The MATLAB codes about the artificial bee_colony (ABC) algorithm%
%          basis_function_coefficients_a :the coefficients%
%          order: the reserved order%
%          Tdata: the duration%



clear; clc;

% Start timing the algorithm
tic;

% Global variables for basis_function_coefficients ranges, time data, order, and degrees of freedom
global basis_function_coefficients_min basis_function_coefficients_max
global time_data truncation_order num_degrees_freedom

% Initialize degrees of freedom and truncation order
num_degrees_freedom = 2;
truncation_order = 49;
time_data = 0:0.01:200;
x=truncation_order + 1;
% Initialize basis_function_coefficients array
basis_function_coefficients = zeros(truncation_order + 1, 1);
basis_function_coefficients(1, 1) = 1;  % Initial frequency
basis_function_coefficients(2, 1) = -0.1;
basis_function_coefficients(3, 1) = 0.1;
initial_basis_function_coefficients = basis_function_coefficients;

% Define basis_function_coefficients range
basis_function_coefficients_min = -3;
basis_function_coefficients_max = 3;

% Bee algorithm basis_function_coefficientss
Optimalvalue = 0;   % Variable for the optimal value
NP = 200;           % Total number of bees
FoodNumber = NP / 2;  % Number of foraging bees
SearchNumber = 15;    % Number of scout bees
maxCycle = 6000;      % Maximum number of iterations
limit = 10;           % Limit for abandoning a source
CR = 0.8;             % Crossover probability for DE algorithm

% Initialize bee positions
Foods = zeros(NP, truncation_order + 1);
for i = 1:NP
    Foods(i, :) = generate_random_solution;
end

% Initialize fitness and trial counters
Fitness = zeros(NP, 1);
trial = zeros(NP, 1);
Foodspath = Foods;

% Calculate initial fitness
for i = 1:NP
    basis_function_coefficients = convert_to_ten(Foods(i, :));
    Fitness(i) = calculate_fitness(basis_function_coefficients);
end

% Sort fitness and find the optimal value
Indfit = sort(Fitness);
Optimalvalue = Indfit(1);
iter = 1;
Zuiyou = zeros(maxCycle, 1);
time_record = zeros(maxCycle, 1);
time_initial = toc;

% Main optimization loop
while (iter <= maxCycle)
    tic;

    % Foraging bee phase
    for i = 1:(FoodNumber)
        % For each foraging bee, a parameter column is randomly selected for optimization
        param_to_change = fix(rand * (truncation_order + 1)) + 1;

        % Randomly select a neighboring solution for comparison
        neighbour = fix(rand * FoodNumber) + 1;
        while (neighbour == i)
            neighbour = fix(rand * FoodNumber) + 1;
        end

        % Copy the current solution
        solution = Foods(i, :);

        % Modify the selected parameter of the solution based on the neighboring solution
        solution(param_to_change) = Foods(i, param_to_change) + ...
            (Foods(i, param_to_change) - Foods(neighbour, param_to_change)) * (rand - 0.5) * 2;

        % Convert the solution to parameter form and calculate its fitness
        converted_parameters = convert_to_ten(solution);
        fitness_sol = calculate_fitness(converted_parameters);

        % Apply greedy selection to update the solution if it has improved
        if (fitness_sol < Fitness(i))
            Foods(i, :) = solution;
            Fitness(i) = fitness_sol;
            trial(i) = 0;
        else
            trial(i) = trial(i) + 1;
            % If a solution hasn't improved for 'limit' trials, generate a new random solution
            if trial(i) > limit
                Foods(i, :) = generate_random_solution;
            end
        end
    end

    % Tournament selection mechanism for onlooker bees
    prob_fit = zeros(1, NP);
    for i = 1:NP
        for j = 1:NP
            if (Fitness(i) > Fitness(j))
                prob_fit(i) = prob_fit(i) + 1;
            end
        end
    end
    prob = (0.9 * prob_fit / max(prob_fit)) + 0.1;

    % Onlooker bee phase
    t = FoodNumber;
    i = 1;
    while (t < (NP - SearchNumber))
        t = t + 1;
        if (rand > prob(i))
            % Select a parameter to modify
            param_to_change = fix(rand * (truncation_order + 1)) + 1;
            neighbour = FoodNumber + fix(rand * (NP - FoodNumber - SearchNumber)) + 1;
            while (neighbour == t)
                neighbour = FoodNumber + fix(rand * (NP - FoodNumber - SearchNumber)) + 1;
            end

            % Copy the solution and modify the selected parameter
            solution = Foods(t, :);
            solution(param_to_change) = Foods(i, param_to_change) + ...
                (Foods(i, param_to_change) - Foods(neighbour, param_to_change)) * (rand - 0.5) * 2;

            % Calculate fitness and apply greedy selection
            converted_parameters = convert_to_ten(solution);
            fitness_sol = calculate_fitness(converted_parameters);
            if (fitness_sol < Fitness(t))
                Foods(t, :) = solution;
                Fitness(t) = fitness_sol;
                trial(t) = 0;
            else
                trial(t) = trial(t) + 1;
                if trial(t) > limit
                    Foods(t, :) = generate_random_solution;
                end
            end
        end
        i = i + 1;
        if i == (FoodNumber + 1)
            i = 1;
        end
    end

    % Scout bee phase
    for i = (NP - SearchNumber + 1):NP
        Foods(i, :) = generate_random_solution;
        converted_parameters = convert_to_ten(Foods(i, :));
        Fitness(i) = calculate_fitness(converted_parameters);
    end


    % Update the global optimal value
    ind = find(Fitness == min(Fitness));
    Min = Fitness(ind);
    if (Min < Optimalvalue)
        Optimalvalue = Min;
        GlobalParams = Foods(ind, :);
    end

    % Output the iteration number and the optimal value
    fprintf('iteration=%d, Optimalvalue=%d\n', iter, Optimalvalue);

    % Record the optimal value and the time for each iteration
    Zuiyou(iter) = Optimalvalue;
    iter = iter + 1;
    toc;
    time_record(iter) = toc;
    Foodspath = [Foodspath Foods];
end

% Calculate total execution time
timesum = sum(time_record) + time_initial

% Plotting the results
figure;
plot(1:maxCycle, Zuiyou, 'r-', 'LineWidth', 1)
xlabel('Iteration Number')
ylabel('Function Value')

