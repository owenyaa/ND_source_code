clear;clc;close all;
%% different initial values
% 不同系统系需要更改的参数
% load 'N_5.mat';
% subplot(3,1,1);
% plot(Tdata,residua(:,1),'k-','LineWidth',1.5);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
% load 'N_15.mat';
% subplot(3,1,2);
% plot(Tdata,residua(:,1),'k-','LineWidth',1.5);
% h1=legend('$$Residues(t)$$');
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
% load 'N_50.mat';
% subplot(3,1,3);
% plot(Tdata,residua(:,1),'k-','LineWidth',1.5);
% h1=legend('$$Non-dimensional Residues(t)$$');
% set(h1,'Interpreter','latex','FontSize',15);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);


load 'N_5.mat';
subplot(3,1,1);
Tdata=0:0.1:10;
a0=parameter_a(1,1);
x=zeros(1,length(Tdata));dx=zeros(1,length(Tdata));
x(1,:)=a0;
for i=1:order
    x(1,:)=x(1,:)+parameter_a(i+1,1)./(1+Tdata).^i;
    dx(1,:)=dx(1,:)-i*parameter_a(i+1,1)./(1+Tdata).^(i+1);
end
y=tanh(Tdata);
%figure;
plot(Tdata,x(1,:),'k-','LineWidth',1);
hold on;
plot(Tdata,y,'r.','MarkerSize',8);
h1=legend('$$TAM$$','$$tanh(t)$$');
set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
load 'N_15.mat';
subplot(3,1,2);
Tdata=0:0.1:10;
a0=parameter_a(1,1);
x=zeros(1,length(Tdata));dx=zeros(1,length(Tdata));
x(1,:)=a0;
for i=1:order
    x(1,:)=x(1,:)+parameter_a(i+1,1)./(1+Tdata).^i;
    dx(1,:)=dx(1,:)-i*parameter_a(i+1,1)./(1+Tdata).^(i+1);
end
y=tanh(Tdata);
%figure;
plot(Tdata,x(1,:),'k-','LineWidth',1);
hold on;
plot(Tdata,y,'r.','MarkerSize',8);
h1=legend('$$TAM$$','$$tanh(t)$$');
set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
load 'N_50.mat';
subplot(3,1,3);
Tdata=0:0.1:10;
a0=parameter_a(1,1);
x=zeros(1,length(Tdata));dx=zeros(1,length(Tdata));
x(1,:)=a0;
for i=1:order
    x(1,:)=x(1,:)+parameter_a(i+1,1)./(1+Tdata).^i;
    dx(1,:)=dx(1,:)-i*parameter_a(i+1,1)./(1+Tdata).^(i+1);
end
y=tanh(Tdata);
% figure;
plot(Tdata,x(1,:),'k-','LineWidth',1);
hold on;
plot(Tdata,y,'r.','MarkerSize',8);
h1=legend('$$TAM$$','$$tanh(t)$$');
set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);




