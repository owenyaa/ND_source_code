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
%          The basis function in this program is fractional function S(t)=\sum_{i=0}%
%          ^{+\infty}a_i\left(1+t\right)^{-it}. The objective function consists of two parts,%
%          the control equation R1 and initial value condition R2.%
%          parameter_a :the coefficients%
%          order: the reserved order%
%          Tdata: the duration%
%          Etol:the convergence error%
%          SSS: response sensitivity matrix%
%          dR:residual vector%
%          lambda_inverse:regularization parameter%

clear; clc; close all;

% Start timing the execution
tic;

% Declare global variables
global time_data basis_function_coefficients truncation_order num_degrees_freedom

% Initialize the number of degrees of freedom
num_degrees_freedom = 2;

% Set the truncation order of the series
truncation_order = 50;

% Define time data for simulation
time_data = 0:0.01:200;

% Initialize basis function coefficients
basis_function_coefficients = zeros(truncation_order + 1, 1);
basis_function_coefficients(1, 1) = 1;  % Initial frequency
basis_function_coefficients(2, 1) = -0.1;
basis_function_coefficients(3, 1) = 0.1;

% Store initial basis function coefficients
ini_basis_function_coefficients = basis_function_coefficients;

% Define function to check coefficient limits
coefficients_limit_check = @(coefficients)(abs(coefficients) < 10);

% Trust-region algorithm parameters
gamma_trust = 1.414;
rho_trust = 0.5;

% Record values of parameters during iteration
coefficients_record = basis_function_coefficients;
trust_region_record = [];

% Maximum numbers for iterations
max_iterations_rs = 1000;  % Response sensitivity iteration
max_iterations_tr = 20;    % Trust region iteration

% Relative error tolerance for convergence
Etol = 1e-10;

% Main loop for response sensitivity iteration
for iteration_index = 1:max_iterations_rs
    % Compute residual identification
    residual_iden = calculate_residual(basis_function_coefficients);

    % Initialize and compute sensitivity matrix
    sensitivity_matrix = reshape(residual_iden(:, num_degrees_freedom + 1:2 * num_degrees_freedom), num_degrees_freedom * length(time_data), 1);
    for i = 1:truncation_order
        sensitivity_matrix = [sensitivity_matrix, reshape(residual_iden(:, num_degrees_freedom * (i + 1) + 1:num_degrees_freedom * (i + 2)), num_degrees_freedom * length(time_data), 1)];
    end

    % Calculate negative residual
    neg_residual = -[residual_iden(:, 1); residual_iden(:, 2)];

    % Perform singular value decomposition
    [U, s, V] = csvd(sensitivity_matrix);

    % Find optimal lambda using L-curve
    lambda_inverse = l_curve(U, s, neg_residual);

    % Initialize temporary coefficients
    temp_coefficients = basis_function_coefficients;

    % Trust-region algorithm loop
    for trust_index = 1:max_iterations_tr
        % Calculate real coefficient updates using Tikhonov regularization
        real_coeff_updates = tikhonov(U, s, V, neg_residual, lambda_inverse);
        temp_real_freq = real_coeff_updates(end, 1);
        sensitivity_coefficients_update = [temp_real_freq; real_coeff_updates(1:end-1, 1)];

        % Check if updated coefficients are within the acceptable range
        if ~coefficients_limit_check(temp_coefficients + sensitivity_coefficients_update)
            lambda_inverse = lambda_inverse * gamma_trust;
            continue;
        end

        % Update basis function coefficients
        basis_function_coefficients = temp_coefficients + sensitivity_coefficients_update;

        % Calculate residuals with new coefficients
        residual_update = calculate_residual(basis_function_coefficients);
        neg_residual_update = -[residual_update(:, 1); residual_update(:, 2)];

        % Calculate agreement indicator
        L_neg_residual = sensitivity_matrix * real_coeff_updates - neg_residual;
        rho_agreement = (neg_residual' * neg_residual - neg_residual_update' * neg_residual_update) / (neg_residual' * neg_residual - L_neg_residual' * L_neg_residual);
        if rho_agreement >= rho_trust
            break;
        end
        lambda_inverse = lambda_inverse * gamma_trust;
    end

    % Check convergence
    tol_convergence = norm(real_coeff_updates) / norm(basis_function_coefficients)
    coefficients_record = [coefficients_record, basis_function_coefficients];
    trust_region_record = [trust_region_record; lambda_inverse];

    if tol_convergence <= Etol
        break;
    end
end

% End timing the execution
toc;

% Calculate final residual
final_residual = calculate_residual(basis_function_coefficients);

% Plot the residuals over time
figure;
plot(time_data, final_residual(:, 1), 'r-', 'LineWidth', 1);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);

% Initialize arrays for computing the system response
system_response = zeros(1, length(time_data));
system_derivative = zeros(1, length(time_data));

% Compute the system response and its derivative
system_response(1, :) = basis_function_coefficients(1, 1);
for i = 1:truncation_order
    system_response(1, :) = system_response(1, :) + basis_function_coefficients(i + 1, 1) ./ (1 + time_data).^i;
    system_derivative(1, :) = system_derivative(1, :) - i * basis_function_coefficients(i + 1, 1) ./ (1 + time_data).^(i + 1);
end

% Compute and plot the tanh function for comparison
tanh_function = tanh(time_data);
figure;
plot(time_data, system_response(1, :), 'k-', 'LineWidth', 1);
hold on;
plot(time_data, tanh_function, 'r.', 'LineWidth', 1);
legend_handle = legend('$$TAM$$', '$$tanh(t)$$');
set(legend_handle, 'Interpreter', 'latex', 'FontSize', 15);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);

% Store initial basis function coefficients for later reference
ini_basis_function_coefficients = basis_function_coefficients;





