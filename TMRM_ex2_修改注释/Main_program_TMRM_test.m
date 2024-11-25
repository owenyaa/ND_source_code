%% Author : GUANG_LIU  * owenyaa@gmail.com *
% Created Time : 2020-11-12 22:13
% Last Revised : GUANG_LIU ,2020-11-12
% Remark : The structure, function and variables of this program are explained:%
%          This program is to calculate the example 2 in the paper "A New Semi-%
%          analytical technique for Nonlinear Systems Based on Response Sensitivity%
%           Analysis", that is, the semi-analytical solution of equation %
%          \begin{equation}%
%          	\label{eq4.15}%
%          	\begin{cases}%
%          		\begin{aligned}%
%          		  & \ddot{\xi}+0.25\ddot{\alpha}+0.1\dot{\xi}+0.2\xi+0.1Q\alpha=0           \\%
%          		  & 0.25\ddot{\xi}+0.5\ddot{\alpha}+0.1\dot{\alpha}-0.04Q\alpha+f(\alpha)=0 %
%          		\end{aligned}%
%          	\end{cases}%
%          \end{equation}%
%          The structure of the program is as follows:%
%          The basis function in this program is 
%          \begin{equation}
%          	\label{eq4.17}
%          	\begin{cases}
%         		\begin{aligned}
%         		  & x_1\approx x^N_1=\sum_{k=1}^{N}\left[b_k\cos\left((2k-1)\omega t\right)+c_k\sin\left((2k-1)\omega t\right)\right] \\
%         		  & x_2\approx x^N_2=\sum_{k=1}^{N}\left[d_k\cos\left((2k-1)\omega t\right)+e_k\sin\left((2k-1)\omega t\right)\right] 
%         		\end{aligned}
%         	\end{cases}
%          \end{equation}
%          harmonic_coefficients :the coefficients%
%          num_harmonics: the reserved order%
%          Etol:the convergence error%
%          sensitivity_matrix: response sensitivity matrix%
%          neg_residual:residual vector%
%          lambda_inverse:regularization parameter%
clear; clc; close all;

% Start timing the execution
tic;

% Declare global variables
global total_time step_size time_data harmonic_coefficients initial_frequency beta radius num_degrees_freedom num_harmonics

% Initialize degrees of freedom and number of harmonics
num_degrees_freedom = 2;
num_harmonics = 15;

% Set initial frequency
initial_frequency = 0.6;

% Set damping coefficient and radius
beta = 20; 
radius = 0;

% Define time parameters
total_time = 2 * pi / initial_frequency;
step_size = 2 * pi / (initial_frequency * 1000);
time_data = (0:step_size:total_time);

% Initialize harmonic coefficients matrix
harmonic_coefficients = zeros(num_harmonics + 1, 2 * num_degrees_freedom);
harmonic_coefficients(1, 1) = initial_frequency;
harmonic_coefficients(2, :) = [0.3, 0.1, 0.2, 0.1];  % Initial values for the first-order harmonics

% Store initial harmonic coefficients
ini_harmonic_coefficients = harmonic_coefficients;

% Define the number of iterations based on time data length
iteration = length(time_data);

% Function to check if coefficients are within acceptable range
coefficients_check = @(harmonic_coefficients)(abs(harmonic_coefficients) < 0.1);

% Algorithm parameters for trust-region algorithm
gamma_trust = 1.414;
rho_trust = 0.5;

% Record values of parameters during iteration
harmonic_coefficients_record = harmonic_coefficients;
trust_region_record = [];

% Maximum number for response sensitivity iteration
max_iterations_rs = 1000;

% Maximum number for trust region iteration
max_iterations_tr = 20;

% Define the relative error tolerance for convergence of the algorithm
Etol = 1e-14;

% Main loop for response sensitivity iteration
for iteration_index = 1:max_iterations_rs
    % Reload frequency and update time parameters
    initial_frequency = harmonic_coefficients(1, 1);
    total_time = 2 * pi / initial_frequency;
    step_size = 2 * pi / (initial_frequency * 1000);
    time_data = (0:step_size:total_time);

    % Compute residual identification
    residual_identification = calculate_residual(harmonic_coefficients);

    % Reshape and adjust residual sensitivity matrix
    temp_sensitivity_matrix = reshape(residual_identification(:, num_degrees_freedom + 1:2 * num_degrees_freedom), num_degrees_freedom * length(time_data), 1);
    for i = 1:2 * num_harmonics * num_degrees_freedom
        temp_sensitivity_matrix = [temp_sensitivity_matrix, reshape(residual_identification(:, num_degrees_freedom * (i + 1) + 1:num_degrees_freedom * (i + 2)), num_degrees_freedom * length(time_data), 1)];
    end
    sensitivity_matrix(:,1:2)=temp_sensitivity_matrix(:,1:2); sensitivity_matrix(:,3:2*num_harmonics*num_degrees_freedom)=temp_sensitivity_matrix(:,4:end);
    
    % Calculate the negative residual
    neg_residual = -reshape(residual_identification(:, 1:num_degrees_freedom), num_degrees_freedom * length(time_data), 1);

    % Singular value decomposition of the sensitivity matrix
    [U, s, V] = csvd(sensitivity_matrix);

    % Find lambda using L-curve
    lambda_inverse = l_curve(U, s, neg_residual);

    % Initialize temporary harmonic coefficients
    temp_harmonic_coefficients = harmonic_coefficients;

    % Trust-region algorithm loop
    for trust_index = 1:max_iterations_tr
        % Calculate real coefficient updates
        real_coeff_updates = tikhonov(U, s, V, neg_residual, lambda_inverse);
        ini_real_coeff_updates=real_coeff_updates;
        temp_real_frequency = zeros(1, 2 * num_degrees_freedom);
        temp_real_frequency(1, 1) = real_coeff_updates(1, 1);
        real_coeff_updates(1, 1) = real_coeff_updates(2, 1);
        real_coeff_updates(2, 1) = 0;

        % Reshape and adjust real coefficient updates
        coeff_updates = reshape(real_coeff_updates, 2, num_degrees_freedom * num_harmonics);
        coeff_updates = coeff_updates';

        % Construct sensitivity parameter update matrix
        sensitivity_parameter_updates = coeff_updates(1:num_harmonics, 1:2);
        for dof_index = 1:num_degrees_freedom - 1
            sensitivity_parameter_updates = [sensitivity_parameter_updates, coeff_updates(dof_index * num_harmonics + 1:(dof_index + 1) * num_harmonics, 1:2)];
        end
        sensitivity_parameter_updates = [temp_real_frequency; sensitivity_parameter_updates];

        % Check if updated parameters are within the acceptable range
        if ~coefficients_check(temp_harmonic_coefficients + sensitivity_parameter_updates)
            lambda_inverse = lambda_inverse * gamma_trust;
            continue;
        end

        % Calculate real coefficient updates
        %         real_coeff_updates = tikhonov(U, s, V, neg_residual, lambda_inverse);
        %         ini_real_coeff_updates=real_coeff_updates;
        %         temp_real_frequency = zeros(1, 2 * num_degrees_freedom);
        %         temp_real_frequency(1, 1) = real_coeff_updates(1, 1);
        %         real_coeff_updates(1, 1) = real_coeff_updates(2, 1);
        %         real_coeff_updates(2, 1) = 0;
        %
        %         % Reshape and adjust real coefficient updates
        %         coeff_updates = reshape(real_coeff_updates, 2, num_degrees_freedom * num_harmonics);
        %         coeff_updates = coeff_updates';
        %
        %         % Construct sensitivity parameter update matrix
        %         sensitivity_parameter_updates = coeff_updates(1:num_harmonics, 1:2);
        %         for dof_index = 1:num_degrees_freedom - 1
        %             sensitivity_parameter_updates = [sensitivity_parameter_updates, coeff_updates(dof_index * num_harmonics + 1:(dof_index + 1) * num_harmonics, 1:2)];
        %         end
        %         sensitivity_parameter_updates = [temp_real_frequency; sensitivity_parameter_updates];

        % Update harmonic coefficients and recalculate residuals
        harmonic_coefficients = temp_harmonic_coefficients + sensitivity_parameter_updates;
        initial_frequency = harmonic_coefficients(1, 1);
        total_time = 2 * pi / initial_frequency;
        step_size = 2 * pi / (initial_frequency * 1000);
        time_data = (0:step_size:total_time);
        residual_updates = calculate_residual(harmonic_coefficients);

        % Calculate agreement indicator and check for convergence
        neg_residual_temp = -reshape(residual_updates(:, 1:num_degrees_freedom), num_degrees_freedom * length(time_data), 1);
        L_neg_residual = sensitivity_matrix * ini_real_coeff_updates - neg_residual;
        rho_agreement = (neg_residual' * neg_residual - neg_residual_temp' * neg_residual_temp) / (neg_residual' * neg_residual - L_neg_residual' * L_neg_residual);
        if rho_agreement >= rho_trust
            break;
        end
        lambda_inverse = lambda_inverse * gamma_trust;
    end

    % Check for overall convergence
    tol_convergence = norm(coeff_updates) / norm(harmonic_coefficients);
    harmonic_coefficients_record = [harmonic_coefficients_record, harmonic_coefficients];
    trust_region_record = [trust_region_record; lambda_inverse];
    harmonic_coefficients
    if tol_convergence <= Etol
        break;
    end
    all_coefficients(iteration_index).harmonic_coefficients = harmonic_coefficients;
    iteration_index;
end

% End timing the execution
toc;

% Calculate final residual
final_residual = calculate_residual(harmonic_coefficients);

% Plot the residuals over time
figure;
plot(time_data, final_residual(:, 1), 'r-', 'LineWidth', 1);
hold on;
plot(time_data, final_residual(:, 2), 'k-', 'LineWidth', 1);
legend_handle = legend('$$h$$', '$$\alpha$$');
set(legend_handle, 'Interpreter', 'latex', 'FontSize', 15);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);

% Define new time data for response analysis
time_data = 0:0.04:200;
initial_frequency = harmonic_coefficients(1, 1);
harmonic_parameters = harmonic_coefficients(2:end, :);

% Initialize arrays for displacement, velocity, and acceleration
displacement = zeros(num_degrees_freedom, length(time_data));
velocity = zeros(num_degrees_freedom, length(time_data));
acceleration = zeros(num_degrees_freedom, length(time_data));

% Compute displacement, velocity, and acceleration for each degree of freedom
for dof_index = 1:num_degrees_freedom
    for harm_index = 1:num_harmonics
        frequency_term = (2 * harm_index - 1) * initial_frequency * time_data;
        displacement(dof_index, :) = displacement(dof_index, :) + harmonic_parameters(harm_index, 2 * dof_index - 1) * cos(frequency_term) + harmonic_parameters(harm_index, 2 * dof_index) * sin(frequency_term);
        velocity(dof_index, :) = velocity(dof_index, :) - initial_frequency * (2 * harm_index - 1) * harmonic_parameters(harm_index, 2 * dof_index - 1) * sin(frequency_term) + initial_frequency * (2 * harm_index - 1) * harmonic_parameters(harm_index, 2 * dof_index) * cos(frequency_term);
        acceleration(dof_index, :) = acceleration(dof_index, :) - (initial_frequency * (2 * harm_index - 1))^2 * harmonic_parameters(harm_index, 2 * dof_index - 1) * cos(frequency_term) - (initial_frequency * (2 * harm_index - 1))^2 * harmonic_parameters(harm_index, 2 * dof_index) * sin(frequency_term);
    end
end

% Plot displacement response over time
figure;
plot(time_data, displacement(1, :), 'k-', 'LineWidth', 1);
hold on;
plot(time_data, displacement(2, :), 'b-', 'LineWidth', 1);
legend_handle = legend('$$h$$', '$$\alpha$$');
set(legend_handle, 'Interpreter', 'latex', 'FontSize', 15);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);

% Plot phase diagrams for displacement and velocity
figure;
plot(displacement(1, :), velocity(1, :), 'k-', 'LineWidth', 1);
hold on;
plot(displacement(2, :), velocity(2, :), 'b-', 'LineWidth', 1);
legend_handle = legend('$$h$$', '$$\alpha$$');
set(legend_handle, 'Interpreter', 'latex', 'FontSize', 15);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);

% Store initial harmonic coefficients for later reference
ini_harmonic_coefficients 