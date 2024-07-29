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
% ��������
% tf=2*pi;h=0.001;Tdata=(0:h:tf/w0)';
total_time=2*pi/w0;step_size=2*pi/(w0*1000);time_data=(0:0.01:10);
C0_1=[0.342172594555505,-0.000133005831228,0.000115963368886]';S0_1=[0.1,-0.000712982230221,-0.000108784126763]';
C0_2=[0.122467788388526,0.001563371451967,-0.000646804589735]';S0_2=[0.072813577124383,0.011674778196277,0.000571448020811]';
parameter_a=zeros(num_harmonics+1,2*num_degrees_freedom);
parameter_a(1,1)=w0;
parameter_a(2:end,:)=[C0_1,S0_1,C0_2,S0_2];   % 6��2��
%% ����в�
residual=calculate_residual_new(parameter_a);
%% ���Ʋв�����
% figure;
% plot(Tdata,residual(:,1),'r-','LineWidth',1);
% hold on;
% plot(Tdata,residual(:,2),'k-','LineWidth',1);
% h1=legend('$$h$$','$$\alpha$$');
% set(h1,'Interpreter','latex','FontSize',15);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);

%% �˴�Ϊ��������֤���򣬱���Ϊ�ò�ִ��������ȣ�������֤�������Ƿ����
% ����������֮�󣬽�ĳ��������ȥһ��С�������²�������һ��
% �������ý�������ٳ���С����Ϊ�����ȣ���֤��ֵ������Ⱥ�ֱ�Ӽ���������������Ƿ��غ�
% ����Ҫע�⣬��ֽ���϶��ǶԵ�(�����Ĳ������)�����ܳ������ԭ��������ݣ�x_cal������
% ������ڶ�Ӧ����ֽ��(parameter_a1)�е�λ��(x1(1,:))���ٶȣ����ٶȶ�Ӧ��ԭ����(parameter_a)�е�����������x_cal(2,:)
% 
for i=1:2*num_harmonics*num_degrees_freedom
    ddt=0.000001;
    %��ʼ������
    parameter_a=zeros(num_harmonics+1,2*num_degrees_freedom);
    parameter_a(1,1)=w0;
    parameter_a(2:end,:)=[C0_1,S0_1,C0_2,S0_2];   % 6��2��
    
    %% ����ϵ������1������˴������Ĳ��
    sensitivity_parameter_a1=zeros(2*num_harmonics*num_degrees_freedom,1);
    sensitivity_parameter_a1(i,1)=1;
    sensitivity_parameter_a1=reshape(sensitivity_parameter_a1,2,num_harmonics*num_degrees_freedom);sensitivity_parameter_a1=sensitivity_parameter_a1';
    sensitivity_parameter_a=sensitivity_parameter_a1(1:num_harmonics,1:2);
    for num_dof=1:num_degrees_freedom-1
        sensitivity_parameter_a=[sensitivity_parameter_a,sensitivity_parameter_a1(num_dof*num_harmonics+1:(num_dof+1)*num_harmonics,1:2)];%da(N_harm+1:2*N_harm,1:2);
    end
    temp_real_w0=zeros(1,2*num_degrees_freedom);
    sensitivity_parameter_a=[temp_real_w0;sensitivity_parameter_a];
    parameter_a=parameter_a+ddt*sensitivity_parameter_a;% ��������ȥ��
    
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
    parameter_a=[C0_1,S0_1,C0_2,S0_2];   % 6��2��
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
















