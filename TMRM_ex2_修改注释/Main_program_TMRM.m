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
%          parameter_a :the coefficients%
%          N_harm: the reserved order%
%          Tdata: the duration%
%          Etol:the convergence error%
%          SSS: response sensitivity matrix%
%          dR:residual vector%
%          lambda_inverse:regularization parameter%

clear;clc;close all;
tic;
global tf h Tdata parameter_a beta r N_dof N_harm
%% different initial values
% case 1
N_dof=2;N_harm=15;
w0=0.6;
% w0=1.4;
beta=20;r=0;
% 其他参数
tf=2*pi/w0;h=2*pi/(w0*1000);Tdata=(0:h:tf);
%% 第一行存储频率，后面存储谐波系数,每个自由度两列
% parameter_a=[w0,  0,    0,   0;
%             C_11,S_11,C_21,S_21;
%             C_12,S_12,C_22,S_22;
%             C_13,S_13,C_23,S_23;...];
parameter_a=zeros(N_harm+1,2*N_dof);
parameter_a(1,1)=w0;
parameter_a(2,:)=[0.3,0.1,0.2,0.1];% 一阶谐波初值
% parameter_a(2,:)=[-0.5,0.1,2,1];% 一阶谐波初值
ini_parameter_a=parameter_a;
iteration=length(Tdata);
% to indicate where a is in the parametric space
parameter_a_judge=@(parameter_a)(abs(parameter_a)<1);
% load simple_fre_data.mat; % load observed data --- Tdata and Xdata
gammaT=1.414;rhob=0.5; % parameter for trust-region algorithm
parameter_a_record=parameter_a; % To record the values of parameters during iteration, and for each iteration, a is recorded in a single row of a_record
TR_record=[];  % recording the parameters during trust region
%% 本算例采用速度和加速度响应识别，且三个自由度噪声等级不同，线加速度1%，角加速度2%和5%
%% Response sensitivity iteration
Nmax=1000;   % maximum number for response sensitivity iteration
Ntr=20;      % maximum number for trust region iteration
%% response sensitivity Solution by ode45
% NT=length(residual(:,1));
for iii=1:Nmax
    %% 重新载入w0,此处为一个周期内均等选取1K个点
    w0=parameter_a(1,1);
    tf=2*pi/w0;h=2*pi/(w0*1000);Tdata=(0:h:tf);
    % compute response and response sensitivity for each incremental
    Etol=1e-10;  % Relative error tolerance for convergence of the RS algorithm
    %%
    residual_iden=cal_residual(parameter_a);
    %% SSS为位移响应灵敏度矩阵，第一列和第二列为残差响应，频率灵敏度从第三列开始
    % 计算的灵敏度矩阵中包含S_11，但是这里要去掉
%     temp_SSS=[residual_iden(:,3);residual_iden(:,4)];%w0
    temp_SSS=reshape(residual_iden(:,N_dof+1:2*N_dof),N_dof*length(Tdata),1);%w0
    for i=1:2*N_harm*N_dof
        temp_SSS=[temp_SSS,reshape(residual_iden(:,N_dof*(i+1)+1:N_dof*(i+2)),N_dof*length(Tdata),1)];
    end
    SSS(:,1:2)=temp_SSS(:,1:2);SSS(:,3:2*N_harm*N_dof)=temp_SSS(:,4:end);
    
    dR=-reshape(residual_iden(:,1:N_dof),N_dof*length(Tdata),1);
    [U,s,V]=csvd(SSS);
    lambda_inverse=l_curve(U,s,dR);
    atemp=parameter_a;
    % trust-region algorithm
    for trust=1:Ntr
        %% 计算的da排序，为w0,C_11,S_11(缺,不参与迭代),C_12,S_12,按照参数矩阵的形式排列da
        real_da=tikhonov(U,s,V,dR,lambda_inverse);
        temp_real_w0=zeros(1,2*N_dof);% 参数矩阵的第一行
        temp_real_w0(1,1)=real_da(1,1);real_da(1,1)=real_da(2,1);real_da(2,1)=0;
        da=reshape(real_da,2,N_dof*N_harm);da=da';
        sensitivity_parameter_da=da(1:N_harm,1:2);
        for num_dof=1:N_dof-1
            sensitivity_parameter_da=[sensitivity_parameter_da,da(num_dof*N_harm+1:(num_dof+1)*N_harm,1:2)];%da(N_harm+1:2*N_harm,1:2);
        end
        sensitivity_parameter_da=[temp_real_w0;sensitivity_parameter_da];
        %%
        if ~parameter_a_judge(atemp+sensitivity_parameter_da)      % if updated a is not out of the parametric space, then, lambda should be increased until updated a is in ...
            lambda_inverse=lambda_inverse*gammaT; %  update of lambda
            continue;
        end
        %% 计算的da排序，为w0,C_11,S_11(缺,不参与迭代),C_12,S_12,按照参数矩阵的形式排列da
        real_da=tikhonov(U,s,V,dR,lambda_inverse);
        ini_da=real_da;
        temp_real_w0=zeros(1,2*N_dof);% 参数矩阵的第一行
        temp_real_w0(1,1)=real_da(1,1);real_da(1,1)=real_da(2,1);real_da(2,1)=0;
        da=reshape(real_da,2,N_dof*N_harm);da=da';
        sensitivity_parameter_da=da(1:N_harm,1:2);
        for num_dof=1:N_dof-1
            sensitivity_parameter_da=[sensitivity_parameter_da,da(num_dof*N_harm+1:(num_dof+1)*N_harm,1:2)];%da(N_harm+1:2*N_harm,1:2);
        end
        sensitivity_parameter_da=[temp_real_w0;sensitivity_parameter_da];
        %%
        parameter_a=atemp+sensitivity_parameter_da;
        %% 重新计算Tdata
        w0=parameter_a(1,1);
        tf=2*pi/w0;h=2*pi/(w0*1000);Tdata=(0:h:tf);
        %% 用新的parameter_a=atemp+da计算响应
        residual_da=cal_residual(parameter_a);

        dRtemp=-reshape(residual_da(:,1:N_dof),N_dof*length(Tdata),1);
        LdR=SSS*ini_da-dR;
        rhos=(dR'*dR-dRtemp'*dRtemp)/(dR'*dR-LdR'*LdR);  % agreement indicator
        if rhos>=rhob
            break;
        end
        lambda_inverse=lambda_inverse*gammaT;
    end
    tolt=norm(da)/norm(parameter_a)
    parameter_a_record=[parameter_a_record,parameter_a];
    TR_record=[TR_record;lambda_inverse];
    parameter_a;
    if tolt<=Etol
        break;
    end
    every_a(iii).parameter_a=parameter_a;
    iii
end
toc;
residua=cal_residual(parameter_a);

figure;
plot(Tdata,residua(:,1),'r-','LineWidth',1);
hold on;
plot(Tdata,residua(:,2),'k-','LineWidth',1);
h1=legend('$$h$$','$$\alpha$$');
set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);


Tdata=0:0.04:200;
w0=parameter_a(1,1);
Harm_parameter_a=parameter_a(2:end,:);
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
% residual=cal_residual(parameter_a);
% figure;
% plot(Tdata,residual(:,1),'k-','LineWidth',1);
% hold on;
% plot(Tdata,residual(:,2),'b-','LineWidth',1);
% h1=legend('$$h$$','$$\alpha$$');
% set(h1,'Interpreter','latex','FontSize',15);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);

figure;
plot(Tdata,x(1,:),'k-','LineWidth',1);
hold on;
plot(Tdata,x(2,:),'b-','LineWidth',1);
h1=legend('$$h$$','$$\alpha$$');
set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);

figure;
plot(x(1,:),dx(1,:),'k-','LineWidth',1);
hold on;
plot(x(2,:),dx(2,:),'b-','LineWidth',1);
h1=legend('$$h$$','$$\alpha$$');
set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);



ini_parameter_a
