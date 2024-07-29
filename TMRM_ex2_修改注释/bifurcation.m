clear;close all;
global M A beta r
beta=20;r=0;
N_dof=2;
% U=25;
M=[1,0.25;0.25,0.5];
C=[0.1,0;0,0.1];
i=1;

U=4.1:0.01:12;
for l=1:length(U)
    K=[0.2,0.1*U(l);0,0.5-0.04*U(l)];
    A=[zeros(N_dof),eye(N_dof);-M\K,-M\C];
    Tdata=0:0.1:1000;
    % ini_x=[x(1,1);x(2,1);dx(1,1);dx(2,1)];
    ini_x=[-0.1;0;0;0];
    [t,num]=ode45('sub_ode45',Tdata,ini_x);
    
    xmax_2=num(8001:10001,2);
    xmax(i).xmax=getmax(xmax_2);
    xmin(i).xmin=getmin(xmax_2);
    i=i+1;
end
for j=1:i-1
    QQ=U(j)*ones(1,length(xmax(j).xmax));
    plot(QQ,xmax(j).xmax,'b.','MarkerSize',2);
    hold on;
    QQ1=U(j)*ones(1,length(xmin(j).xmin));
    plot(QQ1,xmin(j).xmin,'b.','MarkerSize',2);
    hold on;
end
% h1=legend('$$h$$','$$\alpha$$');
% set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);



