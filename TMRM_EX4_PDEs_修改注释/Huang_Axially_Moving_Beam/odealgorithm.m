%% Author : GUANG_LIU  * owenyaa@gmail.com *
% Created Time : 2016-11-03 17:10
% Last Revised : GUANG_LIU ,2016-11-03
% Remark : 本程序的结构、功能和变量做一下说明%

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

% fs=100;%采样频率
% % 采样频率与时间间隔之间的关系： fs=1/dt
% % 采样定理告诉我们，采样频率要大于信号频率的两倍。 
% N=2^16;  %采样点数2^17
% % N个采样点，经过FFT之后，就可以得到N个点的FFT结果。为了方便进行FFT运算，通常N取2的整数次方。
% % 要精确到xHz，则需要采样长度为1/x秒的信号，并做FFT。
% % 要提高频率分辨率，就需要增加采样点数
% n=0:N-1;
% t=n/fs;  % dt=1/fs 表示时间间隔   fs=1/dt
% y=fft(num(700000:end,1),N);  % 进行fft变换
% % 假设采样频率为Fs，信号频率F，采样点数为N。那么FFT之后结果就是一个为N点的复数。
% % 每一个点就对应着一个频率点。这个点的模值，就是该频率值下的幅度特性。
% % y % 输出y看看fft之后的结果。
% m=abs(y(1:N/2))*2/N; % 求信号的真实幅值
% f=2*pi*n*fs/N;  % 注意 乘2*pi 和不乘2*pi的区别%m=log10(m);
% figure;
% % plot(f(1:N/2)./(2*pi),m(1:N/2),'r-','LineWidth',1);
% plot(f(1:N/2),m(1:N/2),'k-','LineWidth',1);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
% h1=legend('$$Iteration steps$$');
% set(h1,'Interpreter','latex','FontSize',15);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
% % 
% fs=100;%采样频率
% % 采样频率与时间间隔之间的关系： fs=1/dt
% % 采样定理告诉我们，采样频率要大于信号频率的两倍。 
% N=2^16;  %采样点数2^17
% % N个采样点，经过FFT之后，就可以得到N个点的FFT结果。为了方便进行FFT运算，通常N取2的整数次方。
% % 要精确到xHz，则需要采样长度为1/x秒的信号，并做FFT。
% % 要提高频率分辨率，就需要增加采样点数
% n=0:N-1;
% t=n/fs;  % dt=1/fs 表示时间间隔   fs=1/dt
% y=fft(num(700000:end,2),N);  % 进行fft变换
% % 假设采样频率为Fs，信号频率F，采样点数为N。那么FFT之后结果就是一个为N点的复数。
% % 每一个点就对应着一个频率点。这个点的模值，就是该频率值下的幅度特性。
% % y % 输出y看看fft之后的结果。
% m=abs(y(1:N/2))*2/N; % 求信号的真实幅值
% f=2*pi*n*fs/N;  % 注意 乘2*pi 和不乘2*pi的区别%m=log10(m);
% figure;
% % plot(f(1:N/2)./(2*pi),m(1:N/2),'r-','LineWidth',1);
% plot(f(1:N/2),m(1:N/2),'r-','LineWidth',1);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
% h1=legend('$$Iteration steps$$');
% set(h1,'Interpreter','latex','FontSize',15);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);