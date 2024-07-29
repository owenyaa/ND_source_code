%% Author : GUANG_LIU  * owenyaa@gmail.com *
% Created Time : 2020-11-12 22:13
% Last Revised : GUANG_LIU ,2020-11-12
% Remark : The structure, function and variables of this program are explained:%
%          This program is to calculate the example 3 in the paper "A New Semi-%
%          analytical technique for Nonlinear Systems Based on Response Sensitivity%
%           Analysis", that is, the semi-analytical solution of equation %
%          \begin{equation}
%          	\label{eq4.24}
%          	\bm{\bar{M}\ddot{x}+\bar{C}\dot{x}+\bar{K}x+\bar{K}^3x^3}=\bm{\bar{F}(t)}
%          \end{equation}
%          The structure of the program is as follows:%
%          The basis function in this program is 
%          \begin{equation}
%          	\label{eq4.26}
%          	\begin{cases}
%          		\begin{aligned}
%          		  & \dot{x}^N_j=\sum_{k=1}^{N}\left[-a_{jk}(2k-1)\omega \sin\left((2k-1)\omega t\right)+b_{jk}(2k-1)\omega\cos((2k-1)\omega t)\right]                     \\
%          		  & \ddot{x}^N_j=\sum_{k=1}^{N}\left[-a_{jk} \left((2k-1)\omega\right)^2\cos((2k-1)\omega t)-b_{jk}\left((2k-1)\omega\right)^2\sin((2k-1)\omega t)\right] \\
%          		\end{aligned}
%          	\end{cases}
%          \end{equation}
%          parameter_a :the coefficients%
%          N_harm: the reserved order%
%          Tdata: the duration%
%          Etol:the convergence error%
%          SSS: response sensitivity matrix%
%          dR:residual vector%
%          lambda_inverse:regularization parameter%

clear; clc; close all;
tic;
% Define global variables
global total_time step_size damping_ratio_1 damping_ratio_4 num_degrees_freedom num_harmonics time_data
global mass_matrix damping_matrix stiffness_matrix force angular_frequency harmonic_coefficients

% Set parameters
num_degrees_freedom = 5;
num_harmonics = 6;
damping_ratio_1 = 2;
damping_ratio_4 = 5;
force=1.5;
angular_frequency = 0.5;

% Define mass, damping, and stiffness matrices
m1 = 1;m2 = 1;m3 = 1;m4 = 1;m5 = 1;

c1 = 0.02 * pi;c2 = 0.05 * pi;c3 = 0.03 * pi;c4 = 0.07 * pi;
c5 = 0.09 * pi;c6 = 0.01 * pi;c7 = 0.04 * pi;c8 = 0.06 * pi;
c9 = 0.08 * pi;c10 = 0.1 * pi;

k1 = 0.2 * pi^2;k2 = 0.8 * pi^2;k3 = 0.5 * pi^2;
k4 = 0.7 * pi^2;k5 = 1.2 * pi^2;k6 = 0.4 * pi^2;
k7 = 1 * pi^2;k8 = 0.3 * pi^2;k9 = 1.1 * pi^2;
k10 = 0.6 * pi^2;

mass_matrix = [m1, 0, 0, 0, 0;
               0, m2, 0, 0, 0;
               0, 0, m3, 0, 0;
               0, 0, 0, m4, 0;
               0, 0, 0, 0, m5];

damping_matrix = [c2 + c3 + c4, 0, -c4, -c3, 0;
                  0, c7 + c8, -c8, 0, 0;
                  -c4, -c8, c4 + c5 + c6 + c8 + c9 + c10, -c5, -c10;
                  -c3, 0, -c5, c1 + c3 + c5, 0;
                  0, 0, -c10, 0, c10];

stiffness_matrix = [k2 + k3 + k4, 0, -k4, -k3, 0;
                    0, k7 + k8, -k8, 0, 0;
                    -k4, -k8, k4 + k5 + k6 + k8 + k9 + k10, -k5, -k10;
                    -k3, 0, -k5, k1 + k3 + k5, 0;
                    0, 0, -k10, 0, k10];

% Define time parameters
total_time = 2*pi/angular_frequency;
step_size = 2*pi/(angular_frequency*1000);
time_data = (0:step_size:total_time);

% Initialize harmonic coefficients matrix
harmonic_coefficients = zeros(num_harmonics, 2*num_degrees_freedom);

% Set the first row to represent angular frequency
harmonic_coefficients(1,:) = [0, 0.03, 0, 0.01, 0, 0.06, 0, 0.05, 0, 0.07];

% Create a copy of the initial harmonic coefficients matrix for reference
initial_harmonic_coefficients = harmonic_coefficients;

% Define the number of iterations
iterations = length(time_data);

% Define a function to check if parameters are within a specified range
coefficients_check = @(parameter)(abs(parameter) < 2);

% Load observed data (if available)
% load simple_frequency_data.mat; 

% Set trust region algorithm parameters
gamma_trust = 1.414;
rho_trust = 0.5;

% Initialize variables to record parameter values during iterations
harmonic_coefficients_record = harmonic_coefficients;
trust_region_record = [];

% Response sensitivity iteration parameters
max_iterations_rs = 1000;
max_iterations_tr = 20;

% Perform response sensitivity iteration
for iteration_index = 1:max_iterations_rs
    % Compute response and response sensitivity for each incremental
    Etol = 1e-10;
    
    % Calculate residual for the current parameter matrix
    current_residual = calculate_residual(harmonic_coefficients);
    
    % Extract response sensitivity matrix
    sensitivity_parameter_updates = reshape(current_residual(:,num_degrees_freedom+1:2*num_degrees_freedom), num_degrees_freedom*length(time_data), 1);
    
    for i = 1:2*num_harmonics*num_degrees_freedom-1
        sensitivity_parameter_updates = [sensitivity_parameter_updates, reshape(current_residual(:,num_degrees_freedom*(i+1)+1:num_degrees_freedom*(i+2)), num_degrees_freedom*length(time_data), 1)];
    end
    
    % Calculate neg_residual (negative of the first column of residual)
    neg_residual = -reshape(current_residual(:,1:num_degrees_freedom), num_degrees_freedom*length(time_data), 1);
    
    % Compute SVD of sensitivity matrix
    [U, s, V] = csvd(sensitivity_parameter_updates);
    
    % Calculate the optimal lambda using L-curve method
    lambda_optimal = l_curve(U, s, neg_residual);
    
    % Store a copy of the harmonic coefficients matrix for later use
    temp_harmonic_coefficients = harmonic_coefficients;
    
    % Trust-region algorithm
    for trust_index = 1:max_iterations_tr
        % Compute real_coeff_updates (optimal parameter increment)
        real_coeff_updates = tikhonov(U, s, V, neg_residual, lambda_optimal);
        
        % Reshape real_coeff_updates to match the harmonic coefficients matrix structure
        coeff_updates = reshape(real_coeff_updates, 2, num_degrees_freedom*num_harmonics);
        coeff_updates = coeff_updates';
        
        % Extract sensitivity_matrix (formerly sensitivity_parameter_da)
        sensitivity_matrix = coeff_updates(1:num_harmonics, 1:2);
        
        for dof_index = 1:num_degrees_freedom-1
            sensitivity_matrix = [sensitivity_matrix, coeff_updates(dof_index*num_harmonics+1:(dof_index+1)*num_harmonics, 1:2)];
        end
        
        % Check if the updated parameter is within the specified range
        if ~coefficients_check(temp_harmonic_coefficients + sensitivity_matrix)
            % If not, increase lambda until the parameter is within range
            lambda_optimal = lambda_optimal * gamma_trust;
            continue;
        end
        
        % Re-compute real_coeff_updates after adjusting lambda
        real_coeff_updates = tikhonov(U, s, V, neg_residual, lambda_optimal);
        
        % Store the initial coeff_updates
        initial_coeff_updates = real_coeff_updates;
        
        % Reshape real_coeff_updates to match the harmonic coefficients matrix structure
        coeff_updates = reshape(real_coeff_updates, 2, num_degrees_freedom*num_harmonics);
        coeff_updates = coeff_updates';
        
        % Extract sensitivity_matrix (formerly sensitivity_parameter_da)
        sensitivity_matrix = coeff_updates(1:num_harmonics, 1:2);
        
        for dof_index = 1:num_degrees_freedom-1
            sensitivity_matrix = [sensitivity_matrix, coeff_updates(dof_index*num_harmonics+1:(dof_index+1)*num_harmonics, 1:2)];
        end
        
        % Update the harmonic coefficients matrix
        harmonic_coefficients = temp_harmonic_coefficients + sensitivity_matrix;
        
        % Re-calculate residual using the updated harmonic coefficients matrix
        residual_updates = calculate_residual(harmonic_coefficients);
        
        % Calculate neg_residual_temp (negative of the first column of residual_da)
        neg_residual_temp = -reshape(residual_updates(:,1:num_degrees_freedom), num_degrees_freedom*length(time_data), 1);
        
        % Calculate L_neg_residual (dot product of sensitivity matrix and initial_coeff_updates)
        L_neg_residual = sensitivity_parameter_updates * initial_coeff_updates - neg_residual;
        
        % Calculate agreement indicator (rhos)
        rhos = (neg_residual' * neg_residual - neg_residual_temp' * neg_residual_temp) / (neg_residual' * neg_residual - L_neg_residual' * L_neg_residual);
        
        % Check if rhos meets the trust region criteria
        if rhos >= rho_trust
            break;
        end
        
        % Adjust lambda if rhos does not meet the criteria
        lambda_optimal = lambda_optimal * gamma_trust;
    end
    
    % Calculate tolerance for convergence
    tol_convergence = norm(coeff_updates) / norm(harmonic_coefficients)
    
    % Store the harmonic coefficients matrix and lambda_inverse
    harmonic_coefficients_record = [harmonic_coefficients_record, harmonic_coefficients];
    trust_region_record = [trust_region_record; lambda_optimal];
    
    % Display the current harmonic coefficients matrix
    harmonic_coefficients;
    
    % Check for convergence
    if tol_convergence <= Etol
        break;
    end
    
    % Store harmonic coefficients matrix at each iteration
    every_a(iteration_index).harmonic_coefficients = harmonic_coefficients;
    iteration_index
end

toc;

% Calculate the final residual
final_residual = calculate_residual(harmonic_coefficients);

% Plot the residuals
figure;
plot(time_data, final_residual(:,1), 'r-', 'LineWidth', 1.5);
hold on;
plot(time_data, final_residual(:,2), 'b-', 'LineWidth', 1.5);
hold on;
plot(time_data, final_residual(:,3), 'k-', 'LineWidth', 1.5);
hold on;
plot(time_data, final_residual(:,4), 'g-', 'LineWidth', 1.5);
hold on;
plot(time_data, final_residual(:,5), 'p-', 'LineWidth', 1.5);
legend('$$x_1$$', '$$x_2$$', '$$x_3$$', '$$x_4$$', '$$x_5$$', 'Interpreter', 'latex', 'FontSize', 15);
xlabel('Time', 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);
ylabel('Residual', 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);

% Define time points for analysis
time_points_analysis = 0:0.01:150;

% Calculate harmonic parameters
harmonic_parameters = harmonic_coefficients(1:end,:);

% Initialize arrays to store displacement, velocity, and acceleration
displacement = zeros(num_degrees_freedom, length(time_points_analysis));
velocity = zeros(num_degrees_freedom, length(time_points_analysis));
acceleration = zeros(num_degrees_freedom, length(time_points_analysis));

% Compute displacement, velocity, and acceleration
for j = 1:num_degrees_freedom
    for i = 1:num_harmonics
        displacement(j,:) = displacement(j,:) + harmonic_parameters(i,2*j-1) * cos((2*i-1)*angular_frequency*time_points_analysis) ...
                            + harmonic_parameters(i,2*j) * sin((2*i-1)*angular_frequency*time_points_analysis);
        velocity(j,:) = velocity(j,:) - angular_frequency*(2*i-1)*harmonic_parameters(i,2*j-1)*sin((2*i-1)*angular_frequency*time_points_analysis) ...
                        + angular_frequency*(2*i-1)*harmonic_parameters(i,2*j)*cos((2*i-1)*angular_frequency*time_points_analysis);
        acceleration(j,:) = acceleration(j,:) - (angular_frequency*(2*i-1))^2*harmonic_parameters(i,2*j-1)*cos((2*i-1)*angular_frequency*time_points_analysis) ...
                            - (angular_frequency*(2*i-1))^2*harmonic_parameters(i,2*j)*sin((2*i-1)*angular_frequency*time_points_analysis);
    end
end

% Plot displacement vs. time
figure;
plot(time_points_analysis, displacement(1,:), 'r-', 'LineWidth', 1.5);
hold on;
plot(time_points_analysis, displacement(2,:), 'b-', 'LineWidth', 1.5);
hold on;
plot(time_points_analysis, displacement(3,:), 'k-', 'LineWidth', 1.5);
hold on;
plot(time_points_analysis, displacement(4,:), 'g-', 'LineWidth', 1.5);
hold on;
plot(time_points_analysis, displacement(5,:), 'p-', 'LineWidth', 1.5);
h1=legend('$$x_1$$', '$$x_2$$', '$$x_3$$', '$$x_4$$', '$$x_5$$', 'Interpreter', 'latex', 'FontSize', 15);
set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
%xlabel('Time', 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);
%ylabel('Displacement', 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);


initial_harmonic_coefficients

