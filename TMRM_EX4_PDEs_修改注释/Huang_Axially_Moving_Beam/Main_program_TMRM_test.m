%% Author : GUANG_LIU  * owenyaa@gmail.com *
% Created Time : 2020-11-12 22:13
% Last Revised : GUANG_LIU ,2020-11-12
% Remark : The structure, function and variables of this program are explained:%
%          This program is to calculate the example 4 in the paper "A New Semi-%
%          analytical technique for Nonlinear Systems Based on Response Sensitivity%
%           Analysis", that is, the semi-analytical solution of equation %
%           \begin{equation}
%           \label{eq4.36}
%               \bm{\tilde{M}\ddot{q}+\tilde{C}\dot{q}+\tilde{K}q+\tilde{N}(q)}=\bm{\tilde{F}(t)}
%           \end{equation}
%          where $\bm{q}=\left[q_1, q_2\right]^T$, the mass matrix $\bm{\tilde{M}}=\begin{pmatrix}1&0\\0&1\end{pmatrix}$; 
%          the damping matrix $\bm{\tilde{C}}=\begin{pmatrix}\mu_{11}&-\mu_{12}\\\mu_{21}&\mu_{22}\end{pmatrix}$; 
%          the linear stiffness matrix $\bm{\tilde{K}}=\begin{pmatrix}k_{11}&0\\0&k_{21}\end{pmatrix}$; 
%          the nonlinear restoring force  $\bm{\tilde{N}(q)}=\begin{pmatrix}k_{12}q_1q^2_2+k_{13}q^3_1\\k_{22}q_2q^2_1+k_{23}q^3_2\end{pmatrix}$ 
%          and the external force vector $\bm{\tilde{F}(t)}=\begin{pmatrix}f_1\cos(\Omega t)\\f_2\cos(\Omega t)\end{pmatrix}$, 
%          in which $\mu_{12}=\mu_{21}=16 v/3$, $k_{11}=\left(v_{f}^{2} \pi^{2}-v^{2}+1\right) \pi^{2}, \quad k_{21}=4\left(4 v_{f}^{2} \pi^{2}-v^{2}+1\right) \pi^{2}$, 
%          $k_{12}=3 v_{1}^{2} \pi^{4}, \quad k_{13}=k_{12} / 8$,$k_{22}=k_{12}, \quad k_{23}=2 k_{12}$, in which $v_{1}^{2}=1124, v_{f}^{2}=0.03,$ and $v=0.6$. 

%          The structure of the program is as follows:%
%          The basis function in this program is 
%          q_j(t)=\sum_{k=1}^{+\infty} \left[c_{jk}\cos\left(k \Omega t\right)+s_{jk}\sin\left(k \Omega t\right)\right], j=1,2\dots n
%          parameter_a :the coefficients%
%          N_harm: the reserved order%
%          Tdata: the duration%
%          Etol:the convergence error%
%          SSS: response sensitivity matrix%
%          dR:residual vector%
%          lambda_inverse:regularization parameter%

% Clear workspace, command window, and close all figures
clear; clc; close all;

% Start timing the execution
tic;

% Declare global variables
global total_time step_size num_degrees_freedom num_harmonics time_data 
global mass_matrix damping_matrix stiffness_matrix harmonic_coefficients force_amplitude1 force_amplitude2 frequency k12 k13 k22 k23

% Initialize system parameters
num_degrees_freedom = 2;
num_harmonics = 10;
%% the parameter please refer to Huang's paper: page 9, Eq.39 and Eq.40
% "Huang, J.L, Zhu, W.D. A new incremental harmonic balance method with two time scales 
% for quasi-periodic motions of an axially moving beam with internal resonance under single-tone external excitation. 
% Journal of Vibration and Acoustics 139(2):021010 (2017)"
velocity = 0.6;
velocity1_squared = 1124;
velocity_fluctuation_squared = 0.03;
u12 = 16 * velocity / 3; u21 = 16 * velocity / 3;
stiffness11 = (velocity_fluctuation_squared * pi^2 - velocity^2 + 1) * pi^2;
stiffness21 = 4 * (4 * velocity_fluctuation_squared * pi^2 - velocity^2 + 1) * pi^2;
k12 = 3 * velocity1_squared * pi^4; k13 = k12 / 8; k22 = k12; k23 = 2 * k12;
damping11 = 0.04; damping22 = 0.04; force_amplitude1 = 0.0055; force_amplitude2 = 0;
omega1 = 2.82232; frequency = 1.15 * omega1;

% Construct matrices
mass_matrix = [1, 0; 0, 1];
damping_matrix = [damping11, -u12; u21, damping22];
stiffness_matrix = [stiffness11, 0; 0, stiffness21];

% Time settings
total_time = 2 * pi / frequency;
step_size = 2 * pi / (frequency * 1000);
time_data = (0:step_size:4 * total_time);

% Initialize harmonic parameters matrix
harmonic_coefficients = zeros(num_harmonics, 2 * num_degrees_freedom);
harmonic_coefficients(1,:) = [-0.002, 0, -0.00001, -0.003];  % Example of LCO_1

% Store initial harmonic parameters
initial_harmonic_coefficients = harmonic_coefficients;

% Prepare for iterative algorithm
iteration_count = length(time_data);
coefficients_check = @(harmonic_coefficients)(abs(harmonic_coefficients) < 1.5);
gamma_trust = 1.414; rho_trust = 0.5; 
harmonic_coefficients_record = harmonic_coefficients; 
trust_region_record = [];

% Begin response sensitivity iteration
max_iterations_rs = 1000;
max_iterations_tr = 20;
error_tolerance = 1e-10;

for iteration_index = 1:max_iterations_rs
    % Calculate residual and sensitivity for current parameters
    residual_current = calculate_residual(harmonic_coefficients);

    % Construct sensitivity matrix
    sensitivity_matrix = reshape(residual_current(:, num_degrees_freedom + 1:2 * num_degrees_freedom), num_degrees_freedom * length(time_data), 1);
    for i = 1:2 * num_harmonics * num_degrees_freedom - 1
        sensitivity_matrix = [sensitivity_matrix, reshape(residual_current(:, num_degrees_freedom * (i + 1) + 1:num_degrees_freedom * (i + 2)), num_degrees_freedom * length(time_data), 1)];
    end

    % Calculate negative residual
    neg_residual = -reshape(residual_current(:, 1:num_degrees_freedom), num_degrees_freedom * length(time_data), 1);

    % Singular value decomposition for sensitivity matrix
    [U, s, V] = csvd(sensitivity_matrix);

    % Find optimal lambda using L-curve
    lambda_optimal = l_curve(U, s, neg_residual);

    % Initialize temporary parameters
    temp_harmonic_coefficients = harmonic_coefficients;

    % Trust-region algorithm
    for trust_iteration = 1:max_iterations_tr
        % Compute updates for harmonic parameters
        real_coeff_updates = tikhonov(U, s, V, neg_residual, lambda_optimal);
        coeff_updates = reshape(real_coeff_updates, 2, num_degrees_freedom * num_harmonics);
        coeff_updates = coeff_updates';
        sensitivity_parameter_updates = coeff_updates(1:num_harmonics, 1:2);
        for dof_index = 1:num_degrees_freedom - 1
            sensitivity_parameter_updates = [sensitivity_parameter_updates, coeff_updates(dof_index * num_harmonics + 1:(dof_index + 1) * num_harmonics, 1:2)];
        end

        % Check if updated parameters are within acceptable range
        if ~coefficients_check(temp_harmonic_coefficients + sensitivity_parameter_updates)
            lambda_optimal = lambda_optimal * gamma_trust;
            continue;
        end

        real_coeff_updates = tikhonov(U, s, V, neg_residual, lambda_optimal);
        coeff_updates = reshape(real_coeff_updates, 2, num_degrees_freedom * num_harmonics);
        coeff_updates = coeff_updates';
        sensitivity_parameter_updates = coeff_updates(1:num_harmonics, 1:2);
        for dof_index = 1:num_degrees_freedom - 1
            sensitivity_parameter_updates = [sensitivity_parameter_updates, coeff_updates(dof_index * num_harmonics + 1:(dof_index + 1) * num_harmonics, 1:2)];
        end

        % Update parameters and calculate new residual
        harmonic_coefficients = temp_harmonic_coefficients + sensitivity_parameter_updates;
        residual_updated = calculate_residual(harmonic_coefficients);

        % Calculate agreement indicator
        neg_residual_temp = -reshape(residual_updated(:, 1:num_degrees_freedom), num_degrees_freedom * length(time_data), 1);
        L_neg_residual = sensitivity_matrix * real_coeff_updates - neg_residual;
        rho_agreement = (neg_residual' * neg_residual - neg_residual_temp' * neg_residual_temp) / (neg_residual' * neg_residual - L_neg_residual' * L_neg_residual);

        % Check for convergence in trust-region algorithm
        if rho_agreement >= rho_trust
            break;
        end
        lambda_optimal = lambda_optimal * gamma_trust;
    end

    % Check for overall convergence
    tol_convergence = norm(coeff_updates) / norm(harmonic_coefficients)
    harmonic_coefficients_record = [harmonic_coefficients_record, harmonic_coefficients];
    trust_region_record = [trust_region_record; lambda_optimal];
    
    if tol_convergence <= error_tolerance
        break;
    end
    all_parameters(iteration_index).harmonic_parameters = harmonic_coefficients;
end

% End timing the execution
toc;

% Calculate final residual
final_residual = calculate_residual(harmonic_coefficients);

% Plot the residuals over time
figure;
plot(time_data, final_residual(:, 1), 'r-', 'LineWidth', 1.5);
hold on;
plot(time_data, final_residual(:, 2), 'b-', 'LineWidth', 1.5);
legend_handle = legend('$$q_1(t)$$', '$$q_2(t)$$');
set(legend_handle, 'Interpreter', 'latex', 'FontSize', 15);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);

% Define new time data for response analysis
time_points_analysis = 0:0.01:50;
harmonic_parameters = harmonic_coefficients(1:end, :);

% Initialize arrays for displacement, velocity, and acceleration
displacement = zeros(num_degrees_freedom, length(time_points_analysis));
velocity = zeros(num_degrees_freedom, length(time_points_analysis));
acceleration = zeros(num_degrees_freedom, length(time_points_analysis));

% Compute displacement, velocity, and acceleration for each degree of freedom
for dof_index = 1:num_degrees_freedom
    for harm_index = 1:num_harmonics
        frequency_term = (2 * harm_index - 1) * frequency * time_points_analysis;
        displacement(dof_index, :) = displacement(dof_index, :) + harmonic_parameters(harm_index, 2 * dof_index - 1) * cos(frequency_term) + harmonic_parameters(harm_index, 2 * dof_index) * sin(frequency_term);
        velocity(dof_index, :) = velocity(dof_index, :) - frequency * (2 * harm_index - 1) * harmonic_parameters(harm_index, 2 * dof_index - 1) * sin(frequency_term) + frequency * (2 * harm_index - 1) * harmonic_parameters(harm_index, 2 * dof_index) * cos(frequency_term);
        acceleration(dof_index, :) = acceleration(dof_index, :) - (frequency * (2 * harm_index - 1))^2 * harmonic_parameters(harm_index, 2 * dof_index - 1) * cos(frequency_term) - (frequency * (2 * harm_index - 1))^2 * harmonic_parameters(harm_index, 2 * dof_index) * sin(frequency_term);
    end
end

% Plot displacement response over time
figure;
plot(time_points_analysis, displacement(1, :), 'r-', 'LineWidth', 1.5);
hold on;
plot(time_points_analysis, displacement(2, :), 'b-', 'LineWidth', 1.5);
legend_handle = legend('$$q_1(t)$$', '$$q_2(t)$$');
set(legend_handle, 'Interpreter', 'latex', 'FontSize', 15);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);

% Store initial harmonic coefficients for later reference
initial_harmonic_coefficients
