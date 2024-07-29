%% Author : GUANG_LIU  * owenyaa@gmail.com *
% Created Time : 2020-11-12 21:08
% Last Revised : GUANG_LIU ,2020-11-12
% Remark : The structure, function and variables of this program are explained:%
%          This document is about the verification procedure of sensitivity matrix,%
%          which verifies the sensitivity response with respect to the parameters by difference.%
clc; clear; close all;

% Declare global variables
global truncation_order time_data

% Initialize the truncation order of the series and time data
truncation_order = 5;
time_data = 0:0.01:20;

% Define and initialize the basis function coefficients
basis_function_coefficients = (0:1:truncation_order)';
basis_function_coefficients(1, 1) = 1;

% Calculate the residual using the basis function coefficients
initial_residual = calculate_residual(basis_function_coefficients);

% Sensitivity verification using difference approximation
% The process alters each coefficient slightly and recalculates the residual
% The difference between these residuals divided by the small change approximates the sensitivity
% This comparison aims to validate the sensitivity calculation against the direct calculation
%% 此处为灵敏度验证程序，本质为用差分代替灵敏度，可以验证灵敏度是否算错
% 计算正问题之后，将某个参数减去一个小量，用新参数再算一次
% 两次所得结果做差再除以小量即为灵敏度，验证差分的灵敏度和直接计算的灵敏度曲线是否重合
% 但是要注意，差分结果肯定是对的(即下文差量结果)，可能出错的是原程序的内容，x_cal中内容
% 另外关于对应，差分结果(parameter_a1)中的位移(x1(1,:))，速度，加速度对应到原程序(parameter_a)中的灵敏度内容x_cal(2,:)
for i = 1:truncation_order
    % Small change for sensitivity calculation
    ddt = 0.000001;
    
    % Reset and slightly alter the current coefficient
    altered_coefficients = (0:1:truncation_order)';
    altered_coefficients(1, 1) = 1;
    altered_coefficients(i + 1, 1) = altered_coefficients(i + 1, 1) + ddt;
    
    % Recalculate residual with altered coefficients
    altered_residual = calculate_residual(altered_coefficients);
    
    % Compute the approximate sensitivity
    approx_sensitivity = (altered_residual - initial_residual) / ddt;
    
    % Plot the original and approximated sensitivity for comparison
    figure; 
    plot(time_data, initial_residual(:, 2 * i + 2), 'r-*')
    hold on
    plot(time_data, approx_sensitivity(:, 2), 'k-*')
end

% Special case for the first coefficient
for i = 1
    % Small change for sensitivity calculation
    ddt = 0.0001;
    
    % Alter the first coefficient
    altered_first_coefficients = (0:1:truncation_order)';
    altered_first_coefficients(1, 1) = 1 + ddt;
    
    % Recalculate residual with altered first coefficient
    altered_residual_first = calculate_residual(altered_first_coefficients);
    
    % Compute the approximate sensitivity for the first coefficient
    approx_sensitivity_first = (altered_residual_first - initial_residual) / ddt;
    
    % Plot the original and approximated sensitivity for the first coefficient
    figure; 
    plot(time_data, initial_residual(:, end), 'r-')
    hold on
    plot(time_data, approx_sensitivity_first(:, 1), 'k-')
end
