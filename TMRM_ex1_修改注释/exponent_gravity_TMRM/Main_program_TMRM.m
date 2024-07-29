%% Author : GUANG_LIU  * owenyaa@gmail.com *
% Created Time : 2020-11-12 21:32
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
%          The basis function in this program is exponential function \chi_i(t)=e^{-it}%
%          The objective function consists of two parts,%
%          the control equation R1 and initial value condition R2.%
%          parameter_a :the coefficients%
%          order: the reserved order%
%          Tdata: the duration%
%          Etol:the convergence error%
%          SSS: response sensitivity matrix%
%          dR:residual vector%
%          lambda_inverse:regularization parameter%

clear;
clc;close all;
tic;
global time_data basis_function_coefficients truncation_order num_degrees_freedom 
num_degrees_freedom=2;
%% different initial values
truncation_order=50;time_data=0:0.01:200;
%% 第一行存储频率，后面存储谐波系数,每个自由度两列
basis_function_coefficients=zeros(truncation_order+1,1);
basis_function_coefficients(1,1)=1;basis_function_coefficients(2,1)=-0.1;basis_function_coefficients(3,1)=0.1;
ini_parameter_a=basis_function_coefficients;
% to indicate where a is in the parametric space
parameter_a_judge=@(parameter_a)(abs(parameter_a)<1000);
% load simple_fre_data.mat; % load observed data --- Tdata and Xdata
gammaT=1.414;rhob=0.5; % parameter for trust-region algorithm
parameter_a_record=basis_function_coefficients; % To record the values of parameters during iteration, and for each iteration, a is recorded in a single row of a_record
TR_record=[];  % recording the parameters during trust region
%% 本算例采用速度和加速度响应识别，且三个自由度噪声等级不同，线加速度1%，角加速度2%和5%
%% Response sensitivity iteration
Nmax=1000;   % maximum number for response sensitivity iteration
Ntr=20;      % maximum number for trust region iteration
%% response sensitivity Solution by ode45
% NT=length(residual(:,1));
for iii=1:Nmax
    %% 重新载入w0,此处为一个周期内均等选取1K个点
    % compute response and response sensitivity for each incremental
    Etol=1e-10;  % Relative error tolerance for convergence of the RS algorithm
    %%
    residual_iden=calculate_residual(basis_function_coefficients);
    %% SSS为位移响应灵敏度矩阵，第一列和第二列为残差响应，频率灵敏度从第三列开始
    % 计算的灵敏度矩阵中包含S_11，但是这里要去掉   
    SSS=reshape(residual_iden(:,num_degrees_freedom+1:2*num_degrees_freedom),num_degrees_freedom*length(time_data),1);
    for j=1:truncation_order
        SSS=[SSS,reshape(residual_iden(:,num_degrees_freedom*(j+1)+1:num_degrees_freedom*(j+2)),num_degrees_freedom*length(time_data),1)];
    end
    dR=-[residual_iden(:,1);residual_iden(:,2)];
    [U,s,V]=csvd(SSS);
    lambda_inverse=l_curve(U,s,dR);
    atemp=basis_function_coefficients;
    % trust-region algorithm
    for trust=1:Ntr
        %% 计算的da排序，为w0,C_11,S_11(缺,不参与迭代),C_12,S_12,按照参数矩阵的形式排列da
        real_da=tikhonov(U,s,V,dR,lambda_inverse);
        temp_real_w0=real_da(end,1);
        sensitivity_parameter_da=[temp_real_w0;real_da(1:end-1,1)];
        %%
        if ~parameter_a_judge(atemp+sensitivity_parameter_da)      % if updated a is not out of the parametric space, then, lambda should be increased until updated a is in ...
            lambda_inverse=lambda_inverse*gammaT; %  update of lambda
            continue;
        end
        %% 计算的da排序，为w0,C_11,S_11(缺,不参与迭代),C_12,S_12,按照参数矩阵的形式排列da
        real_da=tikhonov(U,s,V,dR,lambda_inverse);
        ini_da=real_da;%ini_da(3,1)=0;
        temp_real_w0=real_da(end,1);
        sensitivity_parameter_da=[temp_real_w0;real_da(1:end-1,1)];
        %%
        basis_function_coefficients=atemp+sensitivity_parameter_da;
        %% 用新的parameter_a=atemp+da计算响应
        residual_da=calculate_residual(basis_function_coefficients);
%         dRtemp=-residual_da(:,1);
        dRtemp=-[residual_da(:,1);residual_da(:,2)];
        LdR=SSS*ini_da-dR;
        rhos=(dR'*dR-dRtemp'*dRtemp)/(dR'*dR-LdR'*LdR);  % agreement indicator
        if rhos>=rhob
            break;
        end
        lambda_inverse=lambda_inverse*gammaT;
    end
    tolt=norm(real_da)/norm(basis_function_coefficients)
    parameter_a_record=[parameter_a_record,basis_function_coefficients];
    TR_record=[TR_record;lambda_inverse];
    basis_function_coefficients;
    if tolt<=Etol
        break;
    end
    every_a(iii).parameter_a=basis_function_coefficients;
    iii
end
toc;
residua=cal_residual(basis_function_coefficients);

subplot(2,1,2);
plot(time_data,residua(:,1),'k-','LineWidth',1.5);
% hold on;
% plot(Tdata,residua(:,2),'k-','LineWidth',1);
h1=legend('$$h$$','$$\alpha$$');
set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);

subplot(2,1,1);
time_data=0:0.1:8;
a0=basis_function_coefficients(1,1);
x=zeros(1,length(time_data));dx=zeros(1,length(time_data));
x(1,:)=a0;
for i=1:truncation_order
    x(1,:)=x(1,:)+basis_function_coefficients(i+1,1).*exp(-i.*time_data);
    dx(1,:)=dx(1,:)-i*basis_function_coefficients(i+1,1).*exp(-i.*time_data);
end
y=tanh(time_data);
%figure;
plot(time_data,x(1,:),'k-','LineWidth',1);
hold on;
plot(time_data,y,'r.','MarkerSize',8);
h1=legend('$$TAM$$','$$tanh(t)$$');
set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);




ini_parameter_a
