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

function residual=cal_residual(parameter_a)
global r1 r4 N_dof N_harm Tdata 
global M C K f w
%% 计算方程残差
x=zeros(N_dof,length(Tdata));dx=zeros(N_dof,length(Tdata));ddx=zeros(N_dof,length(Tdata));
for j=1:N_dof
    for i=1:N_harm   % i=1,3,5
        x(j,:)=x(j,:)+parameter_a(i,2*j-1)*cos((2*i-1)*w*Tdata)+parameter_a(i,2*j)*sin((2*i-1)*w*Tdata);
        dx(j,:)=dx(j,:)-w*(2*i-1)*parameter_a(i,2*j-1)*sin((2*i-1)*w*Tdata)+w*(2*i-1)*parameter_a(i,2*j)*cos((2*i-1)*w*Tdata);
        ddx(j,:)=ddx(j,:)-(w*(2*i-1))^2*parameter_a(i,2*j-1)*cos((2*i-1)*w*Tdata)-(w*(2*i-1))^2*parameter_a(i,2*j)*sin((2*i-1)*w*Tdata);
    end
end
non_matrix=zeros(N_dof,N_dof);
non_matrix(1,1)=r1;non_matrix(4,4)=r4;
f_matrix=zeros(N_dof,length(Tdata));f_matrix(3,:)=f*sin(w*Tdata);
residual(1:N_dof,:)=M*ddx+C*dx+K*x+non_matrix*x.^3-f_matrix;
%% 计算谐波系数的灵敏度
for i=1:2*N_harm*N_dof
    sensitivity_parameter_a1=zeros(2*N_harm*N_dof,1);
    x_a=zeros(N_dof,length(Tdata));dx_a=zeros(N_dof,length(Tdata));ddx_a=zeros(N_dof,length(Tdata));
    sensitivity_parameter_a1(i,1)=1;
    sensitivity_parameter_a1=reshape(sensitivity_parameter_a1,2,N_harm*N_dof);sensitivity_parameter_a1=sensitivity_parameter_a1';
    sensitivity_parameter_a=sensitivity_parameter_a1(1:N_harm,1:2);
    for num_dof=1:N_dof-1
        sensitivity_parameter_a=[sensitivity_parameter_a,sensitivity_parameter_a1(num_dof*N_harm+1:(num_dof+1)*N_harm,1:2)];%da(N_harm+1:2*N_harm,1:2);
    end
    for k=1:N_dof
        for j=1:N_harm
            x_a(k,:)=x_a(k,:)+sensitivity_parameter_a(j,2*k-1)*cos((2*j-1)*w*Tdata)+sensitivity_parameter_a(j,2*k)*sin((2*j-1)*w*Tdata);
            dx_a(k,:)=dx_a(k,:)-w*(2*j-1)*sensitivity_parameter_a(j,2*k-1)*sin((2*j-1)*w*Tdata)+w*(2*j-1)*sensitivity_parameter_a(j,2*k)*cos((2*j-1)*w*Tdata);
            ddx_a(k,:)=ddx_a(k,:)-(w*(2*j-1))^2*sensitivity_parameter_a(j,2*k-1)*cos((2*j-1)*w*Tdata)-(w*(2*j-1))^2*sensitivity_parameter_a(j,2*k)*sin((2*j-1)*w*Tdata);
        end
    end
    residual(N_dof*i+1:N_dof*(i+1),:)=M*ddx_a+C*dx_a+K*x_a+3*non_matrix*(x.^2.*x_a);
end

residual=residual';
