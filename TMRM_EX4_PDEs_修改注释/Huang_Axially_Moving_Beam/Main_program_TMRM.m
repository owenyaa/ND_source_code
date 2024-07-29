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
%          The basis function in this program is 
%          q_j(t)=\sum_{k=1}^{+\infty} \left[c_{jk}\cos\left(k \Omega t\right)+s_{jk}\sin\left(k \Omega t\right)\right], j=1,2\dots n
%          parameter_a :the coefficients%
%          N_harm: the reserved order%
%          Tdata: the duration%
%          Etol:the convergence error%
%          SSS: response sensitivity matrix%
%          dR:residual vector%
%          lambda_inverse:regularization parameter%


clear;clc;close all;
tic;
global tf h N_dof N_harm Tdata 
global M C K parameter_a f1 f2 w k12 k13 k22 k23
N_dof=2;N_harm=10;
%% the parameter please refer to Huang's paper: page 9, Eq.39 and Eq.40
% "Huang, J.L, Zhu, W.D. A new incremental harmonic balance method with two time scales 
% for quasi-periodic motions of an axially moving beam with internal resonance under single-tone external excitation. 
% Journal of Vibration and Acoustics 139(2):021010 (2017)"
v=0.6;
v1_2=1124;%v1_2=v1^2
vf_2=0.03;%vf_2=v_f^2
u12=16*v/3; u21=16*v/3;
k11=(vf_2*pi^2-v^2+1)*pi^2;k21=4*(4*vf_2*pi^2-v^2+1)*pi^2;
k12=3*v1_2*pi^4;k13=k12/8;k22=k12;k23=2*k12;
u11=0.04;u22=0.04;f1=0.0055;f2=0;w1=2.82232;w=1.15*w1;
M=[1,0;0,1];C=[u11,-u12;u21,u22];K=[k11,0;0,k21];

tf=2*pi/w;h=2*pi/(w*1000);Tdata=(0:h:4*tf);
%% 第一行存储频率，后面存储谐波系数,每个自由度两列
% parameter_a=[w0,  0,    0,   0;
%             C_11,S_11,C_21,S_21;
%             C_12,S_12,C_22,S_22;
%             C_13,S_13,C_23,S_23;...];
parameter_a=zeros(N_harm,2*N_dof);
parameter_a(1,:)=[-0.002,0,-0.00001,-0.003];%LCO_1
% parameter_a(1,:)=[0.01,0.003,-0.004,0.001];%LCO_2

ini_parameter_a=parameter_a;
iteration=length(Tdata);
% to indicate where a is in the parametric space
parameter_a_judge=@(parameter_a)(abs(parameter_a)<1.5);
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
    % compute response and response sensitivity for each incremental
    Etol=1e-10;  % Relative error tolerance for convergence of the RS algorithm
    %%
    residual_iden=cal_residual(parameter_a);
    %% SSS为位移响应灵敏度矩阵，第一列和第二列为残差响应，频率灵敏度从第三列开始
    SSS=reshape(residual_iden(:,N_dof+1:2*N_dof),N_dof*length(Tdata),1);%w0
    for i=1:2*N_harm*N_dof-1
        SSS=[SSS,reshape(residual_iden(:,N_dof*(i+1)+1:N_dof*(i+2)),N_dof*length(Tdata),1)];
    end
    dR=-reshape(residual_iden(:,1:N_dof),N_dof*length(Tdata),1);
    [U,s,V]=csvd(SSS);
    lambda_inverse=l_curve(U,s,dR);
    atemp=parameter_a;
    % trust-region algorithm
    for trust=1:Ntr
        %% 计算的da排序，为w0,C_11,S_11(缺,不参与迭代),C_12,S_12,按照参数矩阵的形式排列da
        real_da=tikhonov(U,s,V,dR,lambda_inverse);
        da=reshape(real_da,2,N_dof*N_harm);da=da';
        sensitivity_parameter_da=da(1:N_harm,1:2);
        for num_dof=1:N_dof-1
            sensitivity_parameter_da=[sensitivity_parameter_da,da(num_dof*N_harm+1:(num_dof+1)*N_harm,1:2)];%da(N_harm+1:2*N_harm,1:2);
        end
        %%
        if ~parameter_a_judge(atemp+sensitivity_parameter_da)      % if updated a is not out of the parametric space, then, lambda should be increased until updated a is in ...
            lambda_inverse=lambda_inverse*gammaT; %  update of lambda
            continue;
        end
        %% 计算的da排序，为w0,C_11,S_11(缺,不参与迭代),C_12,S_12,按照参数矩阵的形式排列da
        real_da=tikhonov(U,s,V,dR,lambda_inverse);
        ini_da=real_da;
        da=reshape(real_da,2,N_dof*N_harm);da=da';
        sensitivity_parameter_da=da(1:N_harm,1:2);
        for num_dof=1:N_dof-1
            sensitivity_parameter_da=[sensitivity_parameter_da,da(num_dof*N_harm+1:(num_dof+1)*N_harm,1:2)];%da(N_harm+1:2*N_harm,1:2);
        end
        %%
        parameter_a=atemp+sensitivity_parameter_da;
        %% 重新计算Tdata
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
residual=cal_residual(parameter_a);

figure;
plot(Tdata,residual(:,1),'r-','LineWidth',1.5);
hold on;
plot(Tdata,residual(:,2),'b-','LineWidth',1.5);
h1=legend('$$q_1(0)$$','$$q_2(0)$$');
set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);

Tdata1=0:0.01:50;
Harm_parameter_a=parameter_a(1:end,:);
%% 计算方程残差
x=zeros(N_dof,length(Tdata1));dx=zeros(N_dof,length(Tdata1));ddx=zeros(N_dof,length(Tdata1));
for j=1:N_dof
    for i=1:N_harm   % i=1,3,5
        x(j,:)=x(j,:)+Harm_parameter_a(i,2*j-1)*cos((2*i-1)*w*Tdata1)+Harm_parameter_a(i,2*j)*sin((2*i-1)*w*Tdata1);
        dx(j,:)=dx(j,:)-w*(2*i-1)*Harm_parameter_a(i,2*j-1)*sin((2*i-1)*w*Tdata1)+w*(2*i-1)*Harm_parameter_a(i,2*j)*cos((2*i-1)*w*Tdata1);
        ddx(j,:)=ddx(j,:)-(w*(2*i-1))^2*Harm_parameter_a(i,2*j-1)*cos((2*i-1)*w*Tdata1)-(w*(2*i-1))^2*Harm_parameter_a(i,2*j)*sin((2*i-1)*w*Tdata1);
    end
end

figure;
plot(Tdata1,x(1,:),'r-','LineWidth',1.5);
hold on;
plot(Tdata1,x(2,:),'b-','LineWidth',1.5);
h1=legend('$$x_1$$','$$x_2$$');
set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
% figure;
% plot(x(1,:),dx(1,:),'r-','LineWidth',1.5);
% hold on;
% plot(x(2,:),dx(2,:),'b-','LineWidth',1.5);
% h1=legend('$$x_1$$','$$x_2$$');
% set(h1,'Interpreter','latex','FontSize',15);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);

ini_parameter_a
