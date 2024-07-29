%% Author : GUANG_LIU  * owenyaa@gmail.com *
% Created Time : 2020-11-12 21:08
% Last Revised : GUANG_LIU ,2020-11-12
% Remark : The structure, function and variables of this program are explained:%
%          This document is about the verification procedure of sensitivity matrix,%
%          which verifies the sensitivity response with respect to the parameters by difference.%

clc; clear; close all;

% Global variables
global num_degrees_freedom num_harmonics time_data 
global mass_matrix damping_matrix stiffness_matrix external_force1 external_force2 angular_frequency 
global nonlinear_stiffness_coefficients 

% Initialization
num_degrees_freedom = 2;
num_harmonics = 3;
velocity = 0.6;
velocity_squared = 1124; % velocity1^2
nonlinear_velocity_squared = 0.03; % nonlinear velocity^2
damping_coefficient12 = 16 * velocity / 3;
damping_coefficient21 = damping_coefficient12;
stiffness_coefficient11 = (nonlinear_velocity_squared * pi^2 - velocity^2 + 1) * pi^2;
stiffness_coefficient21 = 4 * (4 * nonlinear_velocity_squared * pi^2 - velocity^2 + 1) * pi^2;
nonlinear_stiffness_coefficients = [3 * velocity_squared * pi^4, 3 * velocity_squared * pi^4 / 8, 3 * velocity_squared * pi^4, 2 * 3 * velocity_squared * pi^4];
damping_coefficient11 = 0.04;
damping_coefficient22 = 0.04;
external_force1 = 0.0055;
external_force2 = 0;
fundamental_frequency = 2.82232;
angular_frequency = 1.15 * fundamental_frequency;

% Matrices definition
mass_matrix = [1, 0; 0, 1];
damping_matrix = [damping_coefficient11, -damping_coefficient12; damping_coefficient21, damping_coefficient22];
stiffness_matrix = [stiffness_coefficient11, 0; 0, stiffness_coefficient21];

% Time settings
time_data = 2 * pi / angular_frequency;
step_size = 2 * pi / (angular_frequency * 1000);
time_data = (0:step_size:4 * time_data);

% Harmonic coefficients initialization
cosine_coefficients = [0.3, -0.2, 0.1; 0.4, 0.5, -0.6]';
sine_coefficients = [0.1, -0.2, -0.3; 0.6, 0.5, 0.4]';
harmonic_coefficients = zeros(num_harmonics, 2 * num_degrees_freedom);
harmonic_coefficients(1:end, :) = [cosine_coefficients, sine_coefficients]; 

% Compute residual
residual = calculate_residual(harmonic_coefficients);

%% 此处为灵敏度验证程序，本质为用差分代替灵敏度，可以验证灵敏度是否算错
% 计算正问题之后，将某个参数减去一个小量，用新参数再算一次
% 两次所得结果做差再除以小量即为灵敏度，验证差分的灵敏度和直接计算的灵敏度曲线是否重合
% 但是要注意，差分结果肯定是对的(即下文差量结果)，可能出错的是原程序的内容，residual中内容
% 另外关于对应，差分结果(harmonic_coefficients_temp)中的位移(residual_difference(1,:))，速度，加速度对应到原程序(harmonic_coefficients)residual(2,:)

% Sensitivity analysis using finite differences
for i = 1:2 * num_harmonics * num_degrees_freedom
    delta = 0.000000001;
    harmonic_coefficients_temp = harmonic_coefficients;
    
    sensitivity_array=zeros(2*num_harmonics*num_degrees_freedom,1);
    sensitivity_array(i,1)=1;
    sensitivity_array=reshape(sensitivity_array,2,num_harmonics*num_degrees_freedom);
    sensitivity_array=sensitivity_array';
    sensitivity_params=sensitivity_array(1:num_harmonics,1:2);
    for dof_index=1:num_degrees_freedom-1
        sensitivity_params=[sensitivity_params,sensitivity_array(dof_index*num_harmonics+1:(dof_index+1)*num_harmonics,1:2)];%da(N_harm+1:2*N_harm,1:2);
    end
    harmonic_coefficients_temp=harmonic_coefficients_temp+delta*sensitivity_params;% 增量加上去了

    %harmonic_coefficients_temp = harmonic_coefficients_temp + delta * reshape_reshape_sensitivity_parameters(i, num_harmonics, num_degrees_freedom, delta);
    residual_difference = (calculate_residual(harmonic_coefficients_temp) - residual) / delta;

    % Plotting
    figure;
    plot(time_data, residual(:, num_degrees_freedom * i + 1), 'r-');
    hold on;
    plot(time_data, residual(:, num_degrees_freedom * i + 2), 'r-');
    hold on;
    plot(time_data, residual_difference(:, 1), 'k-');
    hold on;
    plot(time_data, residual_difference(:, 2), 'k-');
end

