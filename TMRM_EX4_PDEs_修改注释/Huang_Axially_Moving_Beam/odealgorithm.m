%% Author : GUANG_LIU  * owenyaa@gmail.com *
% Created Time : 2016-11-03 17:10
% Last Revised : GUANG_LIU ,2016-11-03
% Remark : ������Ľṹ�����ܺͱ�����һ��˵��%

clear;clc;%close all;
%% 
% tic;
global A f1 f2 w k12 k13 k22 k23 M
v=0.6;
v1_2=1124;%v1_2=v1^2
vf_2=0.03;%vf_2=v_f^2
u12=16*v/3; u21=16*v/3;
k11=(vf_2*pi^2-v^2+1)*pi^2;k21=4*(4*vf_2*pi^2-v^2+1)*pi^2;
k12=3*v1_2*pi^4;k13=k12/8;k22=k12;k23=2*k12;

u11=0.04;u22=0.04;f1=0.0055;f2=0;w1=2.82232;w=1.15*w1;

M=[1,0;0,1];C=[u11,-u12;u21,u22];K=[k11,0;0,k21];
A=[zeros(2),eye(2);-K/M,-C/M];
% odex=[x(1,1);x(2,1);dx(1,1);dx(2,1)];
% odex=[0.001;0.00008;0;0];
% odex=[0.15;0;0;0];
q10=-0.15:0.001:-0.1;%q10=-0.1+0.001:0.001:-0.05;%q10=-0.05+0.001:0.001:0;
%q10=0+0.001:0.001:0.05;%q10=0.05+0.001:0.001:0.1;%q10=0.1+0.001:0.001:0.15;
q20=-0.05:0.001:0.05;
index=zeros(length(q10),length(q20));
for i=1:length(q10)
    for j=1:length(q20)
        tic;
        tt=0:0.01:600;odex=[q10(i);q20(j);0;0];
        options=odeset('RelTol',1e-10,'AbsTol',1e-10);
        [t,num]=ode45('odehomo',tt,odex,options);
        if max(abs(num(50000:end,1)))>0.008
            index(i,j)=1;
        end
        i
        j
        toc;
    end
end
save('q10_fu_0_15_fu_0_1.mat','-v7.3');



% tt=0:0.01:800;
% options=odeset('RelTol',1e-10,'AbsTol',1e-10);
% [t,num]=ode45('odehomo',tt,odex,options);
% toc;
% figure;
% hold on;
% plot(tt,num(:,1),'k-');
% hold on
% plot(tt,num(:,2),'k-');
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
% figure;
% plot(num(70000:end,1),num(70000:end,3),'k-');
% hold on
% plot(num(70000:end,2),num(70000:end,4),'r-');
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);

% fs=100;%����Ƶ��
% % ����Ƶ����ʱ����֮��Ĺ�ϵ�� fs=1/dt
% % ��������������ǣ�����Ƶ��Ҫ�����ź�Ƶ�ʵ������� 
% N=2^16;  %��������2^17
% % N�������㣬����FFT֮�󣬾Ϳ��Եõ�N�����FFT�����Ϊ�˷������FFT���㣬ͨ��Nȡ2�������η���
% % Ҫ��ȷ��xHz������Ҫ��������Ϊ1/x����źţ�����FFT��
% % Ҫ���Ƶ�ʷֱ��ʣ�����Ҫ���Ӳ�������
% n=0:N-1;
% t=n/fs;  % dt=1/fs ��ʾʱ����   fs=1/dt
% y=fft(num(700000:end,1),N);  % ����fft�任
% % �������Ƶ��ΪFs���ź�Ƶ��F����������ΪN����ôFFT֮��������һ��ΪN��ĸ�����
% % ÿһ����Ͷ�Ӧ��һ��Ƶ�ʵ㡣������ģֵ�����Ǹ�Ƶ��ֵ�µķ������ԡ�
% % y % ���y����fft֮��Ľ����
% m=abs(y(1:N/2))*2/N; % ���źŵ���ʵ��ֵ
% f=2*pi*n*fs/N;  % ע�� ��2*pi �Ͳ���2*pi������%m=log10(m);
% figure;
% % plot(f(1:N/2)./(2*pi),m(1:N/2),'r-','LineWidth',1);
% plot(f(1:N/2),m(1:N/2),'k-','LineWidth',1);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
% h1=legend('$$Iteration steps$$');
% set(h1,'Interpreter','latex','FontSize',15);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
% % 
% fs=100;%����Ƶ��
% % ����Ƶ����ʱ����֮��Ĺ�ϵ�� fs=1/dt
% % ��������������ǣ�����Ƶ��Ҫ�����ź�Ƶ�ʵ������� 
% N=2^16;  %��������2^17
% % N�������㣬����FFT֮�󣬾Ϳ��Եõ�N�����FFT�����Ϊ�˷������FFT���㣬ͨ��Nȡ2�������η���
% % Ҫ��ȷ��xHz������Ҫ��������Ϊ1/x����źţ�����FFT��
% % Ҫ���Ƶ�ʷֱ��ʣ�����Ҫ���Ӳ�������
% n=0:N-1;
% t=n/fs;  % dt=1/fs ��ʾʱ����   fs=1/dt
% y=fft(num(700000:end,2),N);  % ����fft�任
% % �������Ƶ��ΪFs���ź�Ƶ��F����������ΪN����ôFFT֮��������һ��ΪN��ĸ�����
% % ÿһ����Ͷ�Ӧ��һ��Ƶ�ʵ㡣������ģֵ�����Ǹ�Ƶ��ֵ�µķ������ԡ�
% % y % ���y����fft֮��Ľ����
% m=abs(y(1:N/2))*2/N; % ���źŵ���ʵ��ֵ
% f=2*pi*n*fs/N;  % ע�� ��2*pi �Ͳ���2*pi������%m=log10(m);
% figure;
% % plot(f(1:N/2)./(2*pi),m(1:N/2),'r-','LineWidth',1);
% plot(f(1:N/2),m(1:N/2),'r-','LineWidth',1);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
% h1=legend('$$Iteration steps$$');
% set(h1,'Interpreter','latex','FontSize',15);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);