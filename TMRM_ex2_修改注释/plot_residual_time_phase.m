clear;close all;
global Tdata parameter_a N_dof N_harm
%% different initial values
% 不同系统系需要更改的参数
N_dof=2;
load 'N4.mat';
w0=parameter_a(1,1);
Harm_parameter_a=parameter_a(2:end,:);
tf=2*pi/w0;h=2*pi/(w0*1000);Tdata=(0:h:tf);

%% 计算方程残差
x=zeros(N_dof,length(Tdata));dx=zeros(N_dof,length(Tdata));ddx=zeros(N_dof,length(Tdata));
% 有X,Y两个自由度
% for i=1:N_harm   % i=1,3,5

residual=cal_residual(parameter_a);
subplot(3,1,1);
plot(Tdata,residual(:,1),'r-','LineWidth',1.5);
hold on;
plot(Tdata,residual(:,2),'b-','LineWidth',1.5);
% h1=legend('$$h$$','$$\alpha$$');
% set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);

load 'N8.mat';tf=2*pi/w0;h=2*pi/(w0*1000);Tdata=(0:h:tf);
residual=cal_residual(parameter_a);
subplot(3,1,2);
plot(Tdata,residual(:,1),'r-','LineWidth',1.5);
hold on;
plot(Tdata,residual(:,2),'b-','LineWidth',1.5);
% h1=legend('$$h$$','$$\alpha$$');
% set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);

load 'N12.mat';tf=2*pi/w0;h=2*pi/(w0*1000);Tdata=(0:h:tf);
residual=cal_residual(parameter_a);
subplot(3,1,3);
plot(Tdata,residual(:,1),'r-','LineWidth',1.5);
hold on;
plot(Tdata,residual(:,2),'b-','LineWidth',1.5);
% h1=legend('$$h$$','$$\alpha$$');
% set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
Tdata=0:0.01:100;
w0=parameter_a(1,1);
Harm_parameter_a=parameter_a(2:end,:);
x=zeros(N_dof,length(Tdata));dx=zeros(N_dof,length(Tdata));ddx=zeros(N_dof,length(Tdata));
for j=1:N_dof
    for i=1:N_harm   % i=1,3,5
        x(j,:)=x(j,:)+Harm_parameter_a(i,2*j-1)*cos((2*i-1)*w0*Tdata)+Harm_parameter_a(i,2*j)*sin((2*i-1)*w0*Tdata);
        dx(j,:)=dx(j,:)-w0*(2*i-1)*Harm_parameter_a(i,2*j-1)*sin((2*i-1)*w0*Tdata)+w0*(2*i-1)*Harm_parameter_a(i,2*j)*cos((2*i-1)*w0*Tdata);
        ddx(j,:)=ddx(j,:)-(w0*(2*i-1))^2*Harm_parameter_a(i,2*j-1)*cos((2*i-1)*w0*Tdata)-(w0*(2*i-1))^2*Harm_parameter_a(i,2*j)*sin((2*i-1)*w0*Tdata);
    end
end

global M A beta r
beta=20;r=0;
N_dof=2;
% U=25;
M=[1,0.25;0.25,0.5];
C=[0.1,0;0,0.1];
U=8;

K=[0.2,0.1*U;0,0.5-0.04*U];
A=[zeros(N_dof),eye(N_dof);-M\K,-M\C];
Tdata1=0:0.01:100;
tf=2*pi/w0;h=2*pi/(w0*1000);Tdata1=(0:h:tf);
ini_x=[x(1,1);x(2,1);dx(1,1);dx(2,1)];
% ini_x=[-0.1;0;0;0];
[t,num]=ode45('sub_ode45',Tdata1,ini_x);

% figure;
% plot(Tdata,x(1,:),'r-','LineWidth',1.5);
% hold on;
% plot(Tdata,x(2,:),'b-','LineWidth',1.5);
% hold on 
% plot(Tdata1(1:15:end),num(1:15:end,1),'k.','MarkerSize',8);
% hold on;
% plot(Tdata1(1:30:end),num(1:30:end,2),'k.','MarkerSize',8);
% 
% h1=legend('$$h$$','$$\alpha$$');
% set(h1,'Interpreter','latex','FontSize',15);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
% 
% figure;
% plot(x(1,:),dx(1,:),'r-','LineWidth',1.5);
% hold on;
% plot(x(2,:),dx(2,:),'b-','LineWidth',1.5);
% hold on 
% plot(num(1:15:end,1),num(1:15:end,3),'k.','MarkerSize',8);
% hold on;
% plot(num(1:15:end,2),num(1:15:end,4),'k.','MarkerSize',8);
% 
% h1=legend('$$h$$','$$\alpha$$');
% set(h1,'Interpreter','latex','FontSize',15);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);

data=num(:,1);
N=length(data);
N_fft=2^10;
Y=fft(data,N_fft);
Pyy=2*abs(Y(1:N_fft/2+1))/N_fft;
f=1/h*(0:N_fft/2)/N_fft;
figure;
plot(f,(Pyy(1:(N_fft/2+1))),'r-')
