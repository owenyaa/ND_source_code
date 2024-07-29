% load 'xx1.mat';
clear;%close all;
global M A r1 r4 f w
r1=2;r4=5;f=1.5;
N_dof=5;
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

A=[zeros(N_dof),eye(N_dof);-M\K,-M\C];
for jj=1:1:99
    w=jj*0.1;
    w1(jj)=w;
    Tdata=0:0.01:60*pi/w;
%     ini_x=[x(1,1);x(2,1);x(3,1);x(4,1);x(5,1);dx(1,1);dx(2,1);dx(3,1);dx(4,1);dx(5,1)];
    ini_x=[-0.1;0;0;0;0;0;0;0;0;0];
    [t,num]=ode45('sub_ode45',Tdata,ini_x);
    xmax_1=num(ceil(length(Tdata)*0.9):end,1);xmax_2=num(ceil(length(Tdata)*0.9):end,2);
    xmax_3=num(ceil(length(Tdata)*0.9):end,3);xmax_4=num(ceil(length(Tdata)*0.9):end,4);xmax_5=num(ceil(length(Tdata)*0.9):end,5);
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
end
for i=1:1:99
    temp_xmax_1=xmax1(i).max;
    temp_xmax1(i)=sum(temp_xmax_1)/length(temp_xmax_1);
    temp_xmax_2=xmax2(i).max;
    temp_xmax2(i)=sum(temp_xmax_2)/length(temp_xmax_2);
    temp_xmax_3=xmax3(i).max;
    temp_xmax3(i)=sum(temp_xmax_3)/length(temp_xmax_3);
    temp_xmax_4=xmax4(i).max;
    temp_xmax4(i)=sum(temp_xmax_4)/length(temp_xmax_4);
    temp_xmax_5=xmax5(i).max;
    temp_xmax5(i)=sum(temp_xmax_5)/length(temp_xmax_5);
    temp_xmin_1=xmin1(i).min;
    temp_xmin1(i)=sum(temp_xmin_1)/length(temp_xmin_1);
    temp_xmin_2=xmin2(i).min;
    temp_xmin2(i)=sum(temp_xmin_2)/length(temp_xmin_2);
    temp_xmin_3=xmin3(i).min;
    temp_xmin3(i)=sum(temp_xmin_3)/length(temp_xmin_3);
    temp_xmin_4=xmin4(i).min;
    temp_xmin4(i)=sum(temp_xmin_4)/length(temp_xmin_4);
    temp_xmin_5=xmin5(i).min;
    temp_xmin5(i)=sum(temp_xmin_5)/length(temp_xmin_5);
end
figure;
plot(w1,temp_xmax1,'r*');
hold on;
plot(w1,temp_xmax2,'b*');
hold on;
plot(w1,temp_xmax3,'k*');
hold on;
plot(w1,temp_xmax4,'g*');
hold on;
plot(w1,temp_xmax5,'m*');
figure;
plot(w1,temp_xmin1,'r*');
hold on;
plot(w1,temp_xmin2,'b*');
hold on;
plot(w1,temp_xmin3,'k*');
hold on;
plot(w1,temp_xmin4,'g*');
hold on;
plot(w1,temp_xmin5,'m*');
h1=legend('$$x_1$$','$$x_2$$','$$x_3$$','$$x_4$$','$$x_5$$');
set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);


