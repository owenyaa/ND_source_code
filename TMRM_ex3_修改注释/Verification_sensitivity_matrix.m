%% Author : GUANG_LIU  * owenyaa@gmail.com *
% Created Time : 2020-11-12 21:08
% Last Revised : GUANG_LIU ,2020-11-12
% Remark : The structure, function and variables of this program are explained:%
%          This document is about the verification procedure of sensitivity matrix,%
%          which verifies the sensitivity response with respect to the parameters by difference.%
clc;clear;close all;
global r1 r4 N_dof N_harm Tdata
global M C K f w
N_dof=5;N_harm=3;
r1=2;r4=5;
f=1.5;w=0.5;

m1=1;m2=1;m3=1;m4=1;m5=1;
c1=0.02*pi;c2=0.05*pi;c3=0.03*pi;c4=0.07*pi;c5=0.09*pi;
c6=0.01*pi;c7=0.04*pi;c8=0.06*pi;c9=0.08*pi;c10=0.1*pi;
k1=0.2*pi^2;k2=0.8*pi^2;k3=0.5*pi^2;k4=0.7*pi^2;k5=1.2*pi^2;
k6=0.4*pi^2;k7=1*pi^2;k8=0.3*pi^2;k9=1.1*pi^2;k10=0.6*pi^2;

M=[m1,0,0,0,0;0,m2,0,0,0;0,0,m3,0,0;0,0,0,m4,0;0,0,0,0,m5];
C=[c2+c3+c4,0,-c4,-c3,0;0,c7+c8,-c8,0,0;-c4,-c8,c4+c5+c6+c8+c9+c10,-c5,-c10;...
    -c3,0,-c5,c1+c3+c5,0;0,0,-c10,0,c10];
K=[k2+k3+k4,0,-k4,-k3,0;0,k7+k8,-k8,0,0;-k4,-k8,k4+k5+k6+k8+k9+k10,-k5,-k10;...
    -k3,0,-k5,k1+k3+k5,0;0,0,-k10,0,k10];
% 其他参数
% tf=2*pi;h=0.001;Tdata=(0:h:tf/w0)';
tf=2*pi/w;h=2*pi/(w*1000);Tdata=(0:h:tf);
C0_1=[0.3,-0.2,0.1]';S0_1=[0.1,-0.2,-0.3]';
C0_2=[0.4,0.5,-0.6]';S0_2=[0.6,0.5,0.4]';
C0_3=[0.7,-0.8,0.9]';S0_3=[0.9,-0.8,-0.7]';
C0_4=[0.01,0.02,-0.03]';S0_4=[0.03,0.02,0.01]';
C0_5=[0.04,0.05,-0.06]';S0_5=[0.06,0.05,0.04]';
parameter_a=zeros(N_harm,2*N_dof);
% parameter_a(1,1)=w0;
parameter_a(1:end,:)=[C0_1,S0_1,C0_2,S0_2,C0_3,S0_3,C0_4,S0_4,C0_5,S0_5];   % 6行2列
%% 计算残差
residual=cal_residual(parameter_a);
%% 绘制残差曲线
figure;
plot(Tdata,residual(:,1),'r-','LineWidth',1);
hold on;
plot(Tdata,residual(:,2),'k-','LineWidth',1);
hold on;
plot(Tdata,residual(:,3),'g-','LineWidth',1);
hold on;
plot(Tdata,residual(:,4),'b-','LineWidth',1);
hold on;
plot(Tdata,residual(:,5),'b--','LineWidth',1);
% h1=legend('$$h$$','$$\alpha$$');
% set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);

%% 此处为灵敏度验证程序，本质为用差分代替灵敏度，可以验证灵敏度是否算错
% 计算正问题之后，将某个参数减去一个小量，用新参数再算一次
% 两次所得结果做差再除以小量即为灵敏度，验证差分的灵敏度和直接计算的灵敏度曲线是否重合
% 但是要注意，差分结果肯定是对的(即下文差量结果)，可能出错的是原程序的内容，x_cal中内容
% 另外关于对应，差分结果(parameter_a1)中的位移(x1(1,:))，速度，加速度对应到原程序(parameter_a)中的灵敏度内容x_cal(2,:)
% 
% for i=1:2*N_harm*N_dof
for i=1
    ddt=1;
    %初始化参数
    parameter_a=zeros(N_harm,2*N_dof);
    parameter_a(1:end,:)=[C0_1,S0_1,C0_2,S0_2,C0_3,S0_3,C0_4,S0_4,C0_5,S0_5];   % 6行2列
    
    %% 参数系数矩阵，1代表求此处参数的差分
    sensitivity_parameter_a1=zeros(2*N_harm*N_dof,1);
    sensitivity_parameter_a1(i,1)=1;
    sensitivity_parameter_a1=reshape(sensitivity_parameter_a1,2,N_harm*N_dof);sensitivity_parameter_a1=sensitivity_parameter_a1';
    sensitivity_parameter_a=sensitivity_parameter_a1(1:N_harm,1:2);
    for num_dof=1:N_dof-1
        sensitivity_parameter_a=[sensitivity_parameter_a,sensitivity_parameter_a1(num_dof*N_harm+1:(num_dof+1)*N_harm,1:2)];%da(N_harm+1:2*N_harm,1:2);
    end
    parameter_a=parameter_a+ddt*sensitivity_parameter_a;% 增量加上去了
    
    residual1=cal_residual(parameter_a);
    aaaa=(residual1-residual)/ddt;
    figure; 
    plot(Tdata,residual(:,N_dof*i+1),'r-')
    hold on
    plot(Tdata,residual(:,N_dof*i+2),'r-')
    hold on
    plot(Tdata,residual(:,N_dof*i+3),'r-')
    hold on
    plot(Tdata,residual(:,N_dof*i+4),'r-')
    hold on
    plot(Tdata,residual(:,N_dof*i+5),'r-')
    hold on
    
    plot(Tdata,aaaa(:,1),'k-')
    hold on
    plot(Tdata,aaaa(:,2),'k-')
    hold on
    plot(Tdata,aaaa(:,3),'k-')
    hold on
    plot(Tdata,aaaa(:,4),'k-')
    hold on
    plot(Tdata,aaaa(:,5),'k-')
end


