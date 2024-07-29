%% Author : GUANG_LIU  * owenyaa@gmail.com *
% Created Time : 2020-11-12 21:08
% Last Revised : GUANG_LIU ,2020-11-12
% Remark : The structure, function and variables of this program are explained:%
%          This document is about the verification procedure of sensitivity matrix,%
%          which verifies the sensitivity response with respect to the parameters by difference.%
clc;clear;close all;
global total_time step_size radius beta num_degrees_freedom num_harmonics time_data
%% parameters
% ra=0.5;ah=-0.5;xa=0.25;eta=20;
% gamma=10;
w0=0.681669691477223;
num_degrees_freedom=2;num_harmonics=3;
% real_w0=0.681669691477223;
beta=20;radius=0;
% 其他参数
% tf=2*pi;h=0.001;Tdata=(0:h:tf/w0)';
total_time=2*pi/w0;step_size=2*pi/(w0*1000);time_data=(0:0.01:10);
C0_1=[0.342172594555505,-0.000133005831228,0.000115963368886]';S0_1=[0.1,-0.000712982230221,-0.000108784126763]';
C0_2=[0.122467788388526,0.001563371451967,-0.000646804589735]';S0_2=[0.072813577124383,0.011674778196277,0.000571448020811]';
parameter_a=zeros(num_harmonics+1,2*num_degrees_freedom);
parameter_a(1,1)=w0;
parameter_a(2:end,:)=[C0_1,S0_1,C0_2,S0_2];   % 6行2列
%% 计算残差
residual=calculate_residual_new(parameter_a);
%% 绘制残差曲线
% figure;
% plot(Tdata,residual(:,1),'r-','LineWidth',1);
% hold on;
% plot(Tdata,residual(:,2),'k-','LineWidth',1);
% h1=legend('$$h$$','$$\alpha$$');
% set(h1,'Interpreter','latex','FontSize',15);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);

%% 此处为灵敏度验证程序，本质为用差分代替灵敏度，可以验证灵敏度是否算错
% 计算正问题之后，将某个参数减去一个小量，用新参数再算一次
% 两次所得结果做差再除以小量即为灵敏度，验证差分的灵敏度和直接计算的灵敏度曲线是否重合
% 但是要注意，差分结果肯定是对的(即下文差量结果)，可能出错的是原程序的内容，x_cal中内容
% 另外关于对应，差分结果(parameter_a1)中的位移(x1(1,:))，速度，加速度对应到原程序(parameter_a)中的灵敏度内容x_cal(2,:)
% 
for i=1:2*num_harmonics*num_degrees_freedom
    ddt=0.000001;
    %初始化参数
    parameter_a=zeros(num_harmonics+1,2*num_degrees_freedom);
    parameter_a(1,1)=w0;
    parameter_a(2:end,:)=[C0_1,S0_1,C0_2,S0_2];   % 6行2列
    
    %% 参数系数矩阵，1代表求此处参数的差分
    sensitivity_parameter_a1=zeros(2*num_harmonics*num_degrees_freedom,1);
    sensitivity_parameter_a1(i,1)=1;
    sensitivity_parameter_a1=reshape(sensitivity_parameter_a1,2,num_harmonics*num_degrees_freedom);sensitivity_parameter_a1=sensitivity_parameter_a1';
    sensitivity_parameter_a=sensitivity_parameter_a1(1:num_harmonics,1:2);
    for num_dof=1:num_degrees_freedom-1
        sensitivity_parameter_a=[sensitivity_parameter_a,sensitivity_parameter_a1(num_dof*num_harmonics+1:(num_dof+1)*num_harmonics,1:2)];%da(N_harm+1:2*N_harm,1:2);
    end
    temp_real_w0=zeros(1,2*num_degrees_freedom);
    sensitivity_parameter_a=[temp_real_w0;sensitivity_parameter_a];
    parameter_a=parameter_a+ddt*sensitivity_parameter_a;% 增量加上去了
    
    residual1=calculate_residual_new(parameter_a);
    aaaa=(residual1-residual)/ddt;
    figure; 
    plot(time_data,residual(:,(num_degrees_freedom+1)*(i+1)+1),'r-')
    hold on
    plot(time_data,residual(:,(num_degrees_freedom+1)*(i+1)+2),'r-')
    hold on
    plot(time_data,residual(:,(num_degrees_freedom+1)*(i+2)),'r-')
    hold on
    
    plot(time_data,aaaa(:,1),'k-')
    hold on
    plot(time_data,aaaa(:,2),'k-')
    hold on
    plot(time_data,aaaa(:,3),'k-')
end
% 
% 
for i=1
    % i=8;
    ddt=0.0001;
    w_0=[w0+ddt,0,0,0];
    parameter_a=[C0_1,S0_1,C0_2,S0_2];   % 6行2列
    parameter_a=[w_0;parameter_a];
    residual1=calculate_residual_new(parameter_a);
    aaaa=(residual1-residual)/ddt;
    figure; 
    plot(time_data,residual(:,4),'r-')
    hold on
    plot(time_data,residual(:,5),'r-')
    hold on
    plot(time_data,residual(:,6),'r-')
    hold on
    
    plot(time_data,aaaa(:,1),'k-')
    hold on
    plot(time_data,aaaa(:,2),'k-')
    hold on
    plot(time_data,aaaa(:,3),'k-')
end
















