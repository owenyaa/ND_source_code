%% Author : GUANG_LIU  * owenyaa@gmail.com *
% Created Time : 2020-11-12 21:39
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
%          ^{+\infty}a_i\left(1+t\right)^{-it}. The objective function consists of two parts, 
%          the control equation R1 and initial value condition R2. %
%          %
%          The function of this program is to calculate the residuals, including the%
%          residuals of the control equations, the residuals caused by the initial value
%          conditions, and the sensitivity response of the residuals with %
%          respect to the coefficients.%

function residual = calculate_residual(basis_function_coefficients)
    global truncation_order time_data position velocity

    % Initialize position and velocity as zero arrays
    position = zeros(1, length(time_data));
    velocity = zeros(1, length(time_data));

    % Initial condition based on the first parameter
    initial_condition = basis_function_coefficients(1, 1);
    position(1, :) = initial_condition;

    % Compute the position and velocity using the given parameters
    for i = 1:truncation_order
        position(1, :) = position(1, :) + basis_function_coefficients(i + 1, 1) ./ (1 + time_data).^i;
        velocity(1, :) = velocity(1, :) - i * basis_function_coefficients(i + 1, 1) ./ (1 + time_data).^(i + 1);
    end

    % Calculate the residual for the system
    residual(1, :) = velocity(1, :) + position(1, :).^2 - 1;
    residual(2, :) = sum(basis_function_coefficients);

    % Compute sensitivity of each parameter
    sensitivity_position = zeros(truncation_order, length(time_data));
    sensitivity_velocity = zeros(truncation_order, length(time_data));

    for i = 1:truncation_order
        sensitivity_params = zeros(truncation_order, 1);
        sensitivity_params(i + 1, 1) = 1;  % Set sensitivity for the current parameter

        sensitivity_position(i, :) = sensitivity_position(i, :) + sensitivity_params(i + 1, 1) ./ (1 + time_data).^i;
        sensitivity_velocity(i, :) = sensitivity_velocity(i, :) - i * sensitivity_params(i + 1, 1) ./ (1 + time_data).^(i + 1);

        residual(i * 2 + 1, :) = sensitivity_velocity(i, :) + 2 * position(1, :) .* sensitivity_position(i, :);
        residual((i + 1) * 2, :) = 1;
    end

    % Compute sensitivity of the initial condition
    sensitivity_init_position = ones(1, length(time_data));
    sensitivity_init_velocity = zeros(1, length(time_data));

    residual(2 * (truncation_order + 1) + 1, :) = sensitivity_init_velocity + 2 * position(1, :) .* sensitivity_init_position;
    residual(2 * (truncation_order + 2), :) = 1;

    % Transpose the residual to match expected dimensions
    residual = residual';
end
