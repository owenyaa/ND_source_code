clear;close all;
%global Tdata parameter_a N_dof N_harm
%% different initial values
% 不同系统系需要更改的参数
global tf h r1 r4 N_dof N_harm Tdata
global M C K f w parameter_a
N_dof=5;
r1=2;r4=5;f=1.5;w=0.5;
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
tf=2*pi/w;h=2*pi/(w*1000);Tdata=(0:h:tf);
load 'N1.mat';N_harm=1;
residual=cal_residual(parameter_a);
subplot(3,1,1);
plot(Tdata,residual(:,1),'r-','LineWidth',1.5);
hold on;
plot(Tdata,residual(:,2),'b-','LineWidth',1.5);
hold on;
plot(Tdata,residual(:,3),'k-','LineWidth',1.5);
hold on;
plot(Tdata,residual(:,4),'g-','LineWidth',1.5);
hold on;
plot(Tdata,residual(:,5),'p-','LineWidth',1.5);
% h1=legend('$$x_1$$','$$x_2$$','$$x_3$$','$$x_4$$','$$x_5$$');
% set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);

load 'N3.mat';N_harm=3;
residual=cal_residual(parameter_a);
subplot(3,1,2);
plot(Tdata,residual(:,1),'r-','LineWidth',1.5);
hold on;
plot(Tdata,residual(:,2),'b-','LineWidth',1.5);
hold on;
plot(Tdata,residual(:,3),'k-','LineWidth',1.5);
hold on;
plot(Tdata,residual(:,4),'g-','LineWidth',1.5);
hold on;
plot(Tdata,residual(:,5),'p-','LineWidth',1.5);
% h1=legend('$$x_1$$','$$x_2$$','$$x_3$$','$$x_4$$','$$x_5$$');
% set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);

load 'N6.mat';N_harm=6;
residual=cal_residual(parameter_a);
subplot(3,1,3);
plot(Tdata,residual(:,1),'r-','LineWidth',1.5);
hold on;
plot(Tdata,residual(:,2),'b-','LineWidth',1.5);
hold on;
plot(Tdata,residual(:,3),'k-','LineWidth',1.5);
hold on;
plot(Tdata,residual(:,4),'g-','LineWidth',1.5);
hold on;
plot(Tdata,residual(:,5),'p-','LineWidth',1.5);
% h1=legend('$$x_1$$','$$x_2$$','$$x_3$$','$$x_4$$','$$x_5$$');
% set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);

Tdata=0:0.01:100;
w=parameter_a(1,1);
Harm_parameter_a=parameter_a(2:end,:);
x=zeros(N_dof,length(Tdata));dx=zeros(N_dof,length(Tdata));ddx=zeros(N_dof,length(Tdata));
for j=1:N_dof
    for i=1:N_harm   % i=1,3,5
        x(j,:)=x(j,:)+Harm_parameter_a(i,2*j-1)*cos((2*i-1)*w*Tdata)+Harm_parameter_a(i,2*j)*sin((2*i-1)*w*Tdata);
        dx(j,:)=dx(j,:)-w*(2*i-1)*Harm_parameter_a(i,2*j-1)*sin((2*i-1)*w*Tdata)+w*(2*i-1)*Harm_parameter_a(i,2*j)*cos((2*i-1)*w*Tdata);
        ddx(j,:)=ddx(j,:)-(w*(2*i-1))^2*Harm_parameter_a(i,2*j-1)*cos((2*i-1)*w*Tdata)-(w*(2*i-1))^2*Harm_parameter_a(i,2*j)*sin((2*i-1)*w*Tdata);
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
tf=2*pi/w;h=2*pi/(w*1000);Tdata1=(0:h:tf);
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
