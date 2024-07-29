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
%          The function of this program is to calculate the residuals, including the%
%          residuals of the control equations, the residuals caused by the initial value
%          conditions, and the sensitivity response of the residuals with %
%          respect to the coefficients.%

% function residual=calculate_residual(harmonic_coefficients)
% global damping_ratio_1 damping_ratio_4 num_degrees_freedom num_harmonics time_data 
% global mass_matrix damping_matrix stiffness_matrix fource angular_frequency
% 
% %% 计算方程残差
% x=zeros(num_degrees_freedom,length(time_data));dx=zeros(num_degrees_freedom,length(time_data));ddx=zeros(num_degrees_freedom,length(time_data));
% for j=1:num_degrees_freedom
%     for i=1:num_harmonics   % i=1,3,5
%         x(j,:)=x(j,:)+harmonic_coefficients(i,2*j-1)*cos((2*i-1)*angular_frequency*time_data)+harmonic_coefficients(i,2*j)*sin((2*i-1)*angular_frequency*time_data);
%         dx(j,:)=dx(j,:)-angular_frequency*(2*i-1)*harmonic_coefficients(i,2*j-1)*sin((2*i-1)*angular_frequency*time_data)+angular_frequency*(2*i-1)*harmonic_coefficients(i,2*j)*cos((2*i-1)*angular_frequency*time_data);
%         ddx(j,:)=ddx(j,:)-(angular_frequency*(2*i-1))^2*harmonic_coefficients(i,2*j-1)*cos((2*i-1)*angular_frequency*time_data)-(angular_frequency*(2*i-1))^2*harmonic_coefficients(i,2*j)*sin((2*i-1)*angular_frequency*time_data);
%     end
% end
% non_matrix=zeros(num_degrees_freedom,num_degrees_freedom);
% non_matrix(1,1)=damping_ratio_1;non_matrix(4,4)=damping_ratio_4;
% f_matrix=zeros(num_degrees_freedom,length(time_data));f_matrix(3,:)=fource*sin(angular_frequency*time_data);
% residual(1:num_degrees_freedom,:)=mass_matrix*ddx+damping_matrix*dx+stiffness_matrix*x+non_matrix*x.^3-f_matrix;
% %% 计算谐波系数的灵敏度
% for i=1:2*num_harmonics*num_degrees_freedom
%     sensitivity_parameter_a1=zeros(2*num_harmonics*num_degrees_freedom,1);
%     x_a=zeros(num_degrees_freedom,length(time_data));dx_a=zeros(num_degrees_freedom,length(time_data));ddx_a=zeros(num_degrees_freedom,length(time_data));
%     sensitivity_parameter_a1(i,1)=1;
%     sensitivity_parameter_a1=reshape(sensitivity_parameter_a1,2,num_harmonics*num_degrees_freedom);sensitivity_parameter_a1=sensitivity_parameter_a1';
%     sensitivity_parameter_a=sensitivity_parameter_a1(1:num_harmonics,1:2);
%     for num_dof=1:num_degrees_freedom-1
%         sensitivity_parameter_a=[sensitivity_parameter_a,sensitivity_parameter_a1(num_dof*num_harmonics+1:(num_dof+1)*num_harmonics,1:2)];%da(N_harm+1:2*N_harm,1:2);
%     end
%     for k=1:num_degrees_freedom
%         for j=1:num_harmonics
%             x_a(k,:)=x_a(k,:)+sensitivity_parameter_a(j,2*k-1)*cos((2*j-1)*angular_frequency*time_data)+sensitivity_parameter_a(j,2*k)*sin((2*j-1)*angular_frequency*time_data);
%             dx_a(k,:)=dx_a(k,:)-angular_frequency*(2*j-1)*sensitivity_parameter_a(j,2*k-1)*sin((2*j-1)*angular_frequency*time_data)+angular_frequency*(2*j-1)*sensitivity_parameter_a(j,2*k)*cos((2*j-1)*angular_frequency*time_data);
%             ddx_a(k,:)=ddx_a(k,:)-(angular_frequency*(2*j-1))^2*sensitivity_parameter_a(j,2*k-1)*cos((2*j-1)*angular_frequency*time_data)-(angular_frequency*(2*j-1))^2*sensitivity_parameter_a(j,2*k)*sin((2*j-1)*angular_frequency*time_data);
%         end
%     end
%     residual(num_degrees_freedom*i+1:num_degrees_freedom*(i+1),:)=mass_matrix*ddx_a+damping_matrix*dx_a+stiffness_matrix*x_a+3*non_matrix*(x.^2.*x_a);
% end
% 
% residual=residual';

function residual = calculate_residual(harmonic_coefficients)
global damping_ratio_1 damping_ratio_4 num_degrees_freedom num_harmonics time_data
global mass_matrix damping_matrix stiffness_matrix force angular_frequency

%% Calculate the equation residuals
displacement = zeros(num_degrees_freedom, length(time_data));
velocity = zeros(num_degrees_freedom, length(time_data));
acceleration = zeros(num_degrees_freedom, length(time_data));

% Calculate displacement, velocity, and acceleration for each degree of freedom
for dof = 1:num_degrees_freedom
    for harmonic = 1:num_harmonics
        displacement(dof, :) = displacement(dof, :) + harmonic_coefficients(harmonic, 2 * dof - 1) * cos((2 * harmonic - 1) * angular_frequency * time_data) ...
            + harmonic_coefficients(harmonic, 2 * dof) * sin((2 * harmonic - 1) * angular_frequency * time_data);
        velocity(dof, :) = velocity(dof, :) - angular_frequency * (2 * harmonic - 1) * harmonic_coefficients(harmonic, 2 * dof - 1) * sin((2 * harmonic - 1) * angular_frequency * time_data) ...
            + angular_frequency * (2 * harmonic - 1) * harmonic_coefficients(harmonic, 2 * dof) * cos((2 * harmonic - 1) * angular_frequency * time_data);
        acceleration(dof, :) = acceleration(dof, :) - (angular_frequency * (2 * harmonic - 1))^2 * harmonic_coefficients(harmonic, 2 * dof - 1) * cos((2 * harmonic - 1) * angular_frequency * time_data) ...
            - (angular_frequency * (2 * harmonic - 1))^2 * harmonic_coefficients(harmonic, 2 * dof) * sin((2 * harmonic - 1) * angular_frequency * time_data);
    end
end

damping_matrix_nonlinear = zeros(num_degrees_freedom);
damping_matrix_nonlinear(1, 1) = damping_ratio_1;
damping_matrix_nonlinear(4, 4) = damping_ratio_4;

force_matrix = zeros(num_degrees_freedom, length(time_data));
force_matrix(3, :) = force * sin(angular_frequency * time_data);

% Calculate the equation residuals
residual(1:num_degrees_freedom, :) = mass_matrix * acceleration + damping_matrix * velocity + stiffness_matrix * displacement + damping_matrix_nonlinear * displacement.^3 - force_matrix;

%% Calculate the sensitivity of harmonic coefficients
for i = 1:2 * num_harmonics * num_degrees_freedom
    sensitivity_params = zeros(2 * num_harmonics * num_degrees_freedom, 1);
    displacement_sensitivity = zeros(num_degrees_freedom, length(time_data));
    velocity_sensitivity = zeros(num_degrees_freedom, length(time_data));
    acceleration_sensitivity = zeros(num_degrees_freedom, length(time_data));
    
    sensitivity_params(i, 1) = 1;
    sensitivity_params = reshape(sensitivity_params, 2, num_harmonics * num_degrees_freedom);
    sensitivity_params = sensitivity_params';
    harmonic_sensitivity = sensitivity_params(1:num_harmonics, 1:2);
    
    for num_dof = 1:num_degrees_freedom-1
        harmonic_sensitivity = [harmonic_sensitivity, sensitivity_params(num_dof * num_harmonics + 1:(num_dof + 1) * num_harmonics, 1:2)];
    end
    
    for dof = 1:num_degrees_freedom
        for harmonic = 1:num_harmonics
            displacement_sensitivity(dof, :) = displacement_sensitivity(dof, :) + harmonic_sensitivity(harmonic, 2 * dof - 1) * cos((2 * harmonic - 1) * angular_frequency * time_data) ...
                + harmonic_sensitivity(harmonic, 2 * dof) * sin((2 * harmonic - 1) * angular_frequency * time_data);
            velocity_sensitivity(dof, :) = velocity_sensitivity(dof, :) - angular_frequency * (2 * harmonic - 1) * harmonic_sensitivity(harmonic, 2 * dof - 1) * sin((2 * harmonic - 1) * angular_frequency * time_data) ...
                + angular_frequency * (2 * harmonic - 1) * harmonic_sensitivity(harmonic, 2 * dof) * cos((2 * harmonic - 1) * angular_frequency * time_data);
            acceleration_sensitivity(dof, :) = acceleration_sensitivity(dof, :) - (angular_frequency * (2 * harmonic - 1))^2 * harmonic_sensitivity(harmonic, 2 * dof - 1) * cos((2 * harmonic - 1) * angular_frequency * time_data) ...
                - (angular_frequency * (2 * harmonic - 1))^2 * harmonic_sensitivity(harmonic, 2 * dof) * sin((2 * harmonic - 1) * angular_frequency * time_data);
        end
    end
    
    residual(num_degrees_freedom * i + 1:num_degrees_freedom * (i + 1), :) = mass_matrix * acceleration_sensitivity + damping_matrix * velocity_sensitivity + stiffness_matrix * displacement_sensitivity + 3 * damping_matrix_nonlinear * (displacement.^2 .* displacement_sensitivity);
end

residual = residual';
