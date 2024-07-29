clear;clc;%close all;
global Tdata parameter_a order
%% different initial values
% 不同系统系需要更改的参数
order=50;%Tdata=0:0.01:20;
load 'matlab.mat';
Tdata=0:0.1:6;
a0=parameter_a(1,1);
x=zeros(1,length(Tdata));dx=zeros(1,length(Tdata));
x(1,:)=a0;
for i=1:order
    x(1,:)=x(1,:)+parameter_a(i+1,1).*exp(-i.*Tdata);
    dx(1,:)=dx(1,:)-i*parameter_a(i+1,1).*exp(-i.*Tdata);
end
% residual=cal_residual(parameter_a);
% figure;
% plot(Tdata,residual(:,1),'k-','LineWidth',1);
% hold on;
% plot(Tdata,residual(:,2),'b-','LineWidth',1);
% h1=legend('$$x$$','$$y$$');
% set(h1,'Interpreter','latex','FontSize',15);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
y=tanh(Tdata);
figure;
plot(Tdata,x(1,:),'k-','LineWidth',1.5);
hold on;
plot(Tdata,y,'b.','LineWidth',6);
% 
h1=legend('$$TAM$$','$$tanh(t)$$');
set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);

% figure;
% plot(x(1,:),dx(1,:),'k-','LineWidth',1);
% hold on;
% plot(x(2,:),dx(2,:),'b-','LineWidth',1);
% hold on;
% plot(x(1,:),x(2,:),'r-','LineWidth',1);
% hold on;
% plot(dx(1,:),x(2,:),'k-','LineWidth',1);
% h1=legend('$$x$$','$$y$$','$$\beta$$');
% set(h1,'Interpreter','latex','FontSize',15);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
