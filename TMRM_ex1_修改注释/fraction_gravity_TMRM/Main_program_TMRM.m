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
%          The basis function in this program is fractional function S(t)=\sum_{i=0}%
%          ^{+\infty}a_i\left(1+t\right)^{-it}. The objective function consists of two parts,%
%          the control equation R1 and initial value condition R2.%
%          parameter_a :the coefficients%
%          order: the reserved order%
%          Tdata: the duration%
%          Etol:the convergence error%
%          SSS: response sensitivity matrix%
%          dR:residual vector%
%          lambda_inverse:regularization parameter%

clear;clc;%close all;
tic;
global Tdata parameter_a order N_dof
N_dof=2;
%% different initial values
order=50;Tdata=0:0.01:200;
%% 
parameter_a=zeros(order+1,1);
parameter_a(1,1)=1;parameter_a(2,1)=-0.1;parameter_a(3,1)=0.1;
ini_parameter_a=parameter_a;
% to indicate where a is in the parametric space
parameter_a_judge=@(parameter_a)(abs(parameter_a)<10);
% load simple_fre_data.mat; % load observed data --- Tdata and Xdata
gammaT=1.414;rhob=0.5; % parameter for trust-region algorithm
parameter_a_record=parameter_a; % To record the values of parameters during iteration, and for each iteration, a is recorded in a single row of a_record
TR_record=[];  % recording the parameters during trust region
%% �����������ٶȺͼ��ٶ���Ӧʶ�����������ɶ������ȼ���ͬ���߼��ٶ�1%���Ǽ��ٶ�2%��5%
%% Response sensitivity iteration
Nmax=1000;   % maximum number for response sensitivity iteration
Ntr=20;      % maximum number for trust region iteration
%% response sensitivity Solution by ode45
% NT=length(residual(:,1));
for iii=1:Nmax
    %% ��������w0,�˴�Ϊһ�������ھ���ѡȡ1K����
    % compute response and response sensitivity for each incremental
    Etol=1e-10;  % Relative error tolerance for convergence of the RS algorithm
    %%
    residual_iden=cal_residual(parameter_a);
    %% SSSΪλ����Ӧ�����Ⱦ��󣬵�һ�к͵ڶ���Ϊ�в���Ӧ��Ƶ�������ȴӵ����п�ʼ
    % ����������Ⱦ����а���S_11����������Ҫȥ��   
    SSS=reshape(residual_iden(:,N_dof+1:2*N_dof),N_dof*length(Tdata),1);
    for j=1:order
        SSS=[SSS,reshape(residual_iden(:,N_dof*(j+1)+1:N_dof*(j+2)),N_dof*length(Tdata),1)];
    end
    dR=-[residual_iden(:,1);residual_iden(:,2)];
    [U,s,V]=csvd(SSS);
    lambda_inverse=l_curve(U,s,dR);
    atemp=parameter_a;
    % trust-region algorithm
    for trust=1:Ntr
        %% �����da����Ϊw0,C_11,S_11(ȱ,���������),C_12,S_12,���ղ����������ʽ����da
        real_da=tikhonov(U,s,V,dR,lambda_inverse);
        temp_real_w0=real_da(end,1);
        sensitivity_parameter_da=[temp_real_w0;real_da(1:end-1,1)];
        %%
        if ~parameter_a_judge(atemp+sensitivity_parameter_da)      % if updated a is not out of the parametric space, then, lambda should be increased until updated a is in ...
            lambda_inverse=lambda_inverse*gammaT; %  update of lambda
            continue;
        end
        %% �����da����Ϊw0,C_11,S_11(ȱ,���������),C_12,S_12,���ղ����������ʽ����da
        real_da=tikhonov(U,s,V,dR,lambda_inverse);
        ini_da=real_da;%ini_da(3,1)=0;
        temp_real_w0=real_da(end,1);
        sensitivity_parameter_da=[temp_real_w0;real_da(1:end-1,1)];
        %%
        parameter_a=atemp+sensitivity_parameter_da;
        %% ���µ�parameter_a=atemp+da������Ӧ
        residual_da=cal_residual(parameter_a);
%         dRtemp=-residual_da(:,1);
        dRtemp=-[residual_da(:,1);residual_da(:,2)];
        LdR=SSS*ini_da-dR;
        rhos=(dR'*dR-dRtemp'*dRtemp)/(dR'*dR-LdR'*LdR);  % agreement indicator
        if rhos>=rhob
            break;
        end
        lambda_inverse=lambda_inverse*gammaT;
    end
    tolt=norm(real_da)/norm(parameter_a);
    parameter_a_record=[parameter_a_record,parameter_a];
    TR_record=[TR_record;lambda_inverse];
    parameter_a
    if tolt<=Etol
        break;
    end
    every_a(iii).parameter_a=parameter_a;
    iii
end
toc;
% Tdata=0:0.1:10;
residua=cal_residual(parameter_a);

figure;
plot(Tdata,residua(:,1),'r-','LineWidth',1);
% hold on;
% plot(Tdata,residua(:,2),'k-','LineWidth',1);
% h1=legend('$$h$$','$$\alpha$$');
% set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);

ini_parameter_a


a0=parameter_a(1,1);
x=zeros(1,length(Tdata));dx=zeros(1,length(Tdata));
x(1,:)=a0;
for i=1:order
    x(1,:)=x(1,:)+parameter_a(i+1,1)./(1+Tdata).^i;
    dx(1,:)=dx(1,:)-i*parameter_a(i+1,1)./(1+Tdata).^(i+1);
end
y=tanh(Tdata);
figure;
plot(Tdata,x(1,:),'k-','LineWidth',1);
hold on;
plot(Tdata,y,'r.','LineWidth',1);
h1=legend('$$TAM$$','$$tanh(t)$$');
set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);





