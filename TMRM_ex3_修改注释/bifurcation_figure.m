clear;clc;close all;
global tf h r1 r4 N_dof N_harm Tdata
global M C K f w parameter_a
N_dof=5;N_harm=10;
r1=2;r4=5;f=1.5;%w=0.1;
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
load bifurcation2.mat;
for jj=1:1:991
    w=jj*0.01;
    w1(jj)=w;
    parameter_a=bifurcation_p(jj).parameter_a;
    tf=2*pi/w;h=2*pi/(w*1000);Tdata=(0:h:2*tf);
    %     Tdata=0:0.01:150;
    Harm_parameter_a=parameter_a(1:end,:);
    %% º∆À„∑Ω≥Ã≤–≤Ó
    x=zeros(N_dof,length(Tdata));dx=zeros(N_dof,length(Tdata));ddx=zeros(N_dof,length(Tdata));
    for j=1:N_dof
        for i=1:N_harm   % i=1,3,5
            x(j,:)=x(j,:)+Harm_parameter_a(i,2*j-1)*cos((2*i-1)*w*Tdata)+Harm_parameter_a(i,2*j)*sin((2*i-1)*w*Tdata);
            dx(j,:)=dx(j,:)-w*(2*i-1)*Harm_parameter_a(i,2*j-1)*sin((2*i-1)*w*Tdata)+w*(2*i-1)*Harm_parameter_a(i,2*j)*cos((2*i-1)*w*Tdata);
            ddx(j,:)=ddx(j,:)-(w*(2*i-1))^2*Harm_parameter_a(i,2*j-1)*cos((2*i-1)*w*Tdata)-(w*(2*i-1))^2*Harm_parameter_a(i,2*j)*sin((2*i-1)*w*Tdata);
        end
    end
    x=x';
    xmax_1=x(:,1);xmax_2=x(:,2);xmax_3=x(:,3);xmax_4=x(:,4);xmax_5=x(:,5);
    xmax1(jj).max=getmax(xmax_1);
    xmax2(jj).max=getmax(xmax_2);
    xmax3(jj).max=getmax(xmax_3);
    xmax4(jj).max=getmax(xmax_4);
    xmax5(jj).max=getmax(xmax_5);
    xmin1(jj).min=getmin(xmax_1);
    xmin2(jj).min=getmin(xmax_2);
    xmin3(jj).min=getmin(xmax_3);
    xmin4(jj).min=getmin(xmax_4);
    xmin5(jj).min=getmin(xmax_5);
    jj=jj+1
end
clear;clc;close all;
load matlab.mat;
for i=1:1:991
    xmax_1(i)=xmax1(i).max;
    xmax_2(i)=xmax2(i).max;
    xmax_3(i)=xmax3(i).max;
    xmax_4(i)=xmax4(i).max;
    xmax_5(i)=xmax5(i).max;
    xmin_1(i)=xmin1(i).min;
    xmin_2(i)=xmin2(i).min;
    xmin_3(i)=xmin3(i).min;
    xmin_4(i)=xmin4(i).min;
    xmin_5(i)=xmin5(i).min;
end
figure;
plot(w1,xmax_1,'r-','LineWidth',1.5);
hold on;
plot(w1,xmax_2,'b-','LineWidth',1.5);
hold on;
plot(w1,xmax_3,'k-','LineWidth',1.5);
hold on;
plot(w1,xmax_4,'g-','LineWidth',1.5);
hold on;
plot(w1,xmax_5,'m-','LineWidth',1.5);
% h1=legend('$$x_1$$','$$x_2$$','$$x_3$$','$$x_4$$','$$x_5$$');
% set(h1,'Interpreter','latex','FontSize',15);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
% figure;
hold on;
plot(w1,xmin_1,'r-','LineWidth',1.5);
hold on;
plot(w1,xmin_2,'b-','LineWidth',1.5);
hold on;
plot(w1,xmin_3,'k-','LineWidth',1.5);
hold on;
plot(w1,xmin_4,'g-','LineWidth',1.5);
hold on;
plot(w1,xmin_5,'m-','LineWidth',1.5);
h1=legend('$$x_1$$','$$x_2$$','$$x_3$$','$$x_4$$','$$x_5$$');
set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);

figure;
plot(w1(1:25:end),xmax_1(1:25:end),'r*');
hold on;
plot(w1(1:25:end),xmax_2(1:25:end),'b*');
hold on;
plot(w1(1:25:end),xmax_3(1:25:end),'k*');
hold on;
plot(w1(1:25:end),xmax_4(1:25:end),'g*');
hold on;
plot(w1(1:25:end),xmax_5(1:25:end),'m*');
% h1=legend('$$x_1$$','$$x_2$$','$$x_3$$','$$x_4$$','$$x_5$$');
% set(h1,'Interpreter','latex','FontSize',15);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
% figure;
hold on;
plot(w1(1:25:end),xmin_1(1:25:end),'r*');
hold on;
plot(w1(1:25:end),xmin_2(1:25:end),'b*');
hold on;
plot(w1(1:25:end),xmin_3(1:25:end),'k*');
hold on;
plot(w1(1:25:end),xmin_4(1:25:end),'g*');
hold on;
plot(w1(1:25:end),xmin_5(1:25:end),'m*');
h1=legend('$$x_1$$','$$x_2$$','$$x_3$$','$$x_4$$','$$x_5$$');
set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
