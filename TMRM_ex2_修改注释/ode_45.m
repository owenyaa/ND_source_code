% load 'xx1.mat';
clear;close all;
global M A beta r
beta=20;r=0;
N_dof=2;
% U=25;
M=[1,0.25;0.25,0.5];
C=[0.1,0;0,0.1];
% i=1;
U=8;
% for l=1:length(U)
K=[0.2,0.1*U;0,0.5-0.04*U];
A=[zeros(N_dof),eye(N_dof);-M\K,-M\C];
Tdata=0:0.01:200;
% ini_x=[x(1,1);x(2,1);dx(1,1);dx(2,1)];
ini_x=[-0.1;0;0;0];
[t,num]=ode45('sub_ode45',Tdata,ini_x);
    
%     xmax_2=num(8001:10001,1);
%     xmax(i).xmax=getmax(xmax_2);
%     xmin(i).xmin=getmin(xmax_2);
%     i=i+1;
% end
% for j=1:i-1
%     QQ=U(j)*ones(1,length(xmax(j).xmax));
%     plot(QQ,xmax(j).xmax,'k.','MarkerSize',6);
%     hold on;
%     QQ1=U(j)*ones(1,length(xmin(j).xmin));
%     plot(QQ1,xmin(j).xmin,'b.','MarkerSize',6);
%     hold on;
% end
% 
% data=x_cal(:,1);
% N=length(data);
% N_fft=2^14;
% Y=fft(data,N_fft);
% Pyy=2*abs(Y(1:N_fft/2+1))/N_fft;
% f=1/h*(0:N_fft/2)/N_fft;
% figure;
% plot(f,(Pyy(1:(N_fft/2+1))),'k-')
% 
%     QQ=temp(i).Q;QQ=QQ*ones(1,length(xmax));
%     hold on;
% plot(QQ,xmax,'k.','MarkerSize',20);

%
figure;
plot(Tdata,num(:,1),'r-','LineWidth',1.5);
hold on;
plot(Tdata,num(:,2),'k-','LineWidth',1.5);
h1=legend('$$h$$','$$\alpha$$');
set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
% rectangle('Position',[60,-0.4,90,0.8]);
% rectangle('Position',[70,-0.45,90,0.9]);

% figure;
% plot(Tdata(1:10:end),num(1:10:end,1),'r.','LineWidth',1.5);
% hold on;
% plot(Tdata(1:10:end),num(1:10:end,2),'k.','LineWidth',1.5);
% h1=legend('$$h$$','$$\alpha$$');
% set(h1,'Interpreter','latex','FontSize',15);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);

% figure;
% plot(num(1:50:end,1),num(1:50:end,3),'r.','LineWidth',1.5);
% hold on;
% plot(num(1:50:end,2),num(1:50:end,4),'k.','LineWidth',1.5);
% h1=legend('$$h$$','$$\alpha$$');
% set(h1,'Interpreter','latex','FontSize',15);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);

