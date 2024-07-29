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
%          The function of this program is to calculate the residuals, including the%
%          residuals of the control equations, the residuals caused by the initial value
%          conditions, and the sensitivity response of the residuals with %
%          respect to the coefficients.%
function residual=cal_residual(parameter_a)
global r beta N_dof N_harm Tdata
w0=parameter_a(1,1);
Harm_parameter_a=parameter_a(2:end,:);
% Tdata=(0:0.01:10);
%注意，此处gp和gv符号弄反了
% miu=100;damph=0;dampa=0;freratio=0.25;m=2.1647;%b=0.075
% miu=100;damph=0;dampa=0;gp=0;gv=0;m=2.1647;%b=0.075
%U=1.11.21.31.51.82.0
%% 自由度数目，谐波数目，参数个数2*N_harm*N_dof
% N_dof=2;N_harm=3;
U=8;
M=[1,0.25;0.25,0.5];
C=[0.1,0;0,0.1];
K=[0.2,0.1*U;0,0.5-0.04*U];
%% 计算方程残差
x=zeros(N_dof,length(Tdata));dx=zeros(N_dof,length(Tdata));ddx=zeros(N_dof,length(Tdata));
% 有X,Y两个自由度
% for i=1:N_harm   % i=1,3,5
for j=1:N_dof
    for i=1:N_harm   % i=1,3,5
        x(j,:)=x(j,:)+Harm_parameter_a(i,2*j-1)*cos((2*i-1)*w0*Tdata)+Harm_parameter_a(i,2*j)*sin((2*i-1)*w0*Tdata);
        dx(j,:)=dx(j,:)-w0*(2*i-1)*Harm_parameter_a(i,2*j-1)*sin((2*i-1)*w0*Tdata)+w0*(2*i-1)*Harm_parameter_a(i,2*j)*cos((2*i-1)*w0*Tdata);
        ddx(j,:)=ddx(j,:)-(w0*(2*i-1))^2*Harm_parameter_a(i,2*j-1)*cos((2*i-1)*w0*Tdata)-(w0*(2*i-1))^2*Harm_parameter_a(i,2*j)*sin((2*i-1)*w0*Tdata);
    end
end
residual(1:N_dof,:)=M*ddx+C*dx+K*x+[r,0;0,beta]*x.^3;
%% 计算频率的灵敏度
x_w0=zeros(N_dof,length(Tdata));dx_w0=zeros(N_dof,length(Tdata));ddx_w0=zeros(N_dof,length(Tdata));
% y_w0=0;dy_w0=0;ddy_w0=0;
for k=1:N_dof
    for i=1:N_harm   % i=1,3,5
        x_w0(k,:)=x_w0(k,:)-Harm_parameter_a(i,2*k-1)*(2*i-1)*Tdata.*sin((2*i-1)*w0*Tdata)+Harm_parameter_a(i,2*k)*(2*i-1)*Tdata.*cos((2*i-1)*w0*Tdata);
        dx_w0(k,:)=dx_w0(k,:)-((2*i-1)*Harm_parameter_a(i,2*k-1)*sin((2*i-1)*w0*Tdata)+w0*Tdata*(2*i-1)^2*Harm_parameter_a(i,2*k-1).*cos((2*i-1)*w0*Tdata))+...
            ((2*i-1)*Harm_parameter_a(i,2*k)*cos((2*i-1)*w0*Tdata)-w0*Tdata*(2*i-1)^2*Harm_parameter_a(i,2*k).*sin((2*i-1)*w0*Tdata));
        ddx_w0(k,:)=ddx_w0(k,:)-(2*w0*(2*i-1)^2*Harm_parameter_a(i,2*k-1)*cos((2*i-1)*w0*Tdata)-w0^2*Tdata*(2*i-1)^3*Harm_parameter_a(i,2*k-1).*sin((2*i-1)*w0*Tdata))-...
            (2*w0*(2*i-1)^2*Harm_parameter_a(i,2*k)*sin((2*i-1)*w0*Tdata)+w0^2*Tdata*(2*i-1)^3*Harm_parameter_a(i,2*k).*cos((2*i-1)*w0*Tdata));
    end
end
residual(N_dof+1:2*N_dof,:)=M*ddx_w0+C*dx_w0+K*x_w0+3*[r,0;0,beta]*(x.^2.*x_w0);
%% 计算谐波系数的灵敏度
for i=1:2*N_harm*N_dof
    sensitivity_parameter_a1=zeros(2*N_harm*N_dof,1);
    x_a=zeros(N_dof,length(Tdata));dx_a=zeros(N_dof,length(Tdata));ddx_a=zeros(N_dof,length(Tdata));
    %     y_a=0;dy_a=0;ddy_a=0;
    sensitivity_parameter_a1(i,1)=1;
    sensitivity_parameter_a1=reshape(sensitivity_parameter_a1,2,N_harm*N_dof);sensitivity_parameter_a1=sensitivity_parameter_a1';
    sensitivity_parameter_a=sensitivity_parameter_a1(1:N_harm,1:2);
    for num_dof=1:N_dof-1
        sensitivity_parameter_a=[sensitivity_parameter_a,sensitivity_parameter_a1(num_dof*N_harm+1:(num_dof+1)*N_harm,1:2)];%da(N_harm+1:2*N_harm,1:2);
    end
    for k=1:N_dof
        for j=1:N_harm
            x_a(k,:)=x_a(k,:)+sensitivity_parameter_a(j,2*k-1)*cos((2*j-1)*w0*Tdata)+sensitivity_parameter_a(j,2*k)*sin((2*j-1)*w0*Tdata);
            dx_a(k,:)=dx_a(k,:)-w0*(2*j-1)*sensitivity_parameter_a(j,2*k-1)*sin((2*j-1)*w0*Tdata)+w0*(2*j-1)*sensitivity_parameter_a(j,2*k)*cos((2*j-1)*w0*Tdata);
            ddx_a(k,:)=ddx_a(k,:)-(w0*(2*j-1))^2*sensitivity_parameter_a(j,2*k-1)*cos((2*j-1)*w0*Tdata)-(w0*(2*j-1))^2*sensitivity_parameter_a(j,2*k)*sin((2*j-1)*w0*Tdata);
        end
    end
    residual(N_dof*(i+1)+1:N_dof*(i+2),:)=M*ddx_a+C*dx_a+K*x_a+3*[r,0;0,beta]*(x.^2.*x_a);
end

residual=residual';




