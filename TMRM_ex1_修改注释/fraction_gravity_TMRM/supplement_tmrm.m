clear;clc;close all;
tic;
global Tdata parameter_a order N_dof
N_dof=2;
%% different initial values
order=50;Tdata=0:0.01:200;
%% 第一行存储频率，后面存储谐波系数,每个自由度两列
parameter_a=zeros(order+1,1);
parameter_a(1,1)=1;parameter_a(2,1)=-0.1;parameter_a(3,1)=0.1;
ini_parameter_a=parameter_a;
% to indicate where a is in the parametric space
parameter_a_judge=@(parameter_a)(abs(parameter_a)<10);
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
    Etol=1e-6;  % Relative error tolerance for convergence of the RS algorithm
    %%
    residual_iden=cal_residual(parameter_a);
    %% SSS为位移响应灵敏度矩阵，第一列和第二列为残差响应，频率灵敏度从第三列开始
    % 计算的灵敏度矩阵中包含S_11，但是这里要去掉   
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
        parameter_a=atemp+sensitivity_parameter_da;
        %% 用新的parameter_a=atemp+da计算响应
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
t_record=toc
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





