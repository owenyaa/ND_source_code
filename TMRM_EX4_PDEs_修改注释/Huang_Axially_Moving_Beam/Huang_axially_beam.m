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
%          the response is calculated by the 4-th order RK method

clear; clc; close all;
%% the parameter please refer to Huang's paper: page 9, Eq.39 and Eq.40
% "Huang, J.L, Zhu, W.D. A new incremental harmonic balance method with two time scales 
% for quasi-periodic motions of an axially moving beam with internal resonance under single-tone external excitation. 
% Journal of Vibration and Acoustics 139(2):021010 (2017)"

% Global variables declaration
global system_matrix force_amplitude1 force_amplitude2 angular_frequency nonlinear_stiffness_coeffs mass_matrix

% System parameters initialization (refer to Huang's paper)
velocity = 0.6;
velocity_squared = 1124;
nonlinear_velocity_squared = 0.03;
coupling_damping12 = 16 * velocity / 3; 
coupling_damping21 = 16 * velocity / 3;
stiffness11 = (nonlinear_velocity_squared * pi^2 - velocity^2 + 1) * pi^2;
stiffness21 = 4 * (4 * nonlinear_velocity_squared * pi^2 - velocity^2 + 1) * pi^2;
nonlinear_stiffness_coeffs = [3 * velocity_squared * pi^4, 3 * velocity_squared * pi^4 / 8, 3 * velocity_squared * pi^4, 2 * 3 * velocity_squared * pi^4];
damping11 = 0.04; damping22 = 0.04;

force_amplitude1 = 0.0055; force_amplitude2 = 0;
fundamental_frequency = 2.82232; angular_frequency = 1.15 * fundamental_frequency;

% Matrices definition
mass_matrix = [1, 0; 0, 1];
damping_matrix = [damping11, -coupling_damping12; coupling_damping21, damping22];
stiffness_matrix = [stiffness11, 0; 0, stiffness21];

% System matrix for ODE
system_matrix = [zeros(2), eye(2); -inv(mass_matrix) * stiffness_matrix, -inv(mass_matrix) * damping_matrix];

% Initial conditions and time span
%第一组初值，x2做第一种形式的周期运动
initial_conditions = [0.026; -0.039; 0; 0];
%第二组初值，x2做第二种形式的周期运动
% initial_conditions=[0;0;0;0];

time_span = 0:0.01:500;
ode_options = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);

% Solve ODE
[t, solution] = ode45(@ode_axially_beam, time_span, initial_conditions, ode_options);

% Plot results
figure;
plot(t, solution(:, 1), 'r-', 'LineWidth', 1.5);
hold on;
plot(t, solution(:, 2), 'b-', 'LineWidth', 1.5);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);
legend_handle = legend('$$q_1$$', '$$q_2$$');
set(legend_handle, 'Interpreter', 'latex', 'FontSize', 15);

% Additional plots for specified time ranges and states
figure;
plot(t(190000:end), solution(190000:end, 1), 'k-');
hold on;
plot(t(190000:end), solution(190000:end, 2), 'r-');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);

figure;
plot(solution(180000:end, 1), solution(180000:end, 3), 'k-');
hold on;
plot(solution(180000:end, 2), solution(180000:end, 4), 'r-');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);

function y_prime = ode_axially_beam(t, x)
    % Access global variables
    global system_matrix force_amplitude1 force_amplitude2 angular_frequency nonlinear_stiffness_coeffs mass_matrix

    % Nonlinear force vector
    nonlinear_force = [0; 0; mass_matrix \ [force_amplitude1 * cos(angular_frequency * t) - nonlinear_stiffness_coeffs(1) * x(1) * x(2)^2 - nonlinear_stiffness_coeffs(2) * x(1)^3; 
                                             force_amplitude2 * cos(angular_frequency * t) - nonlinear_stiffness_coeffs(3) * x(2) * x(1)^2 - nonlinear_stiffness_coeffs(4) * x(2)^3]];

    % System ODE
    y_prime = system_matrix * x + nonlinear_force;
end


