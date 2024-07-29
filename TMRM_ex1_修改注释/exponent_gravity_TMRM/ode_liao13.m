clear;clc
ini_x=0;
tt=0:0.1:8;
[t,num]=ode45('sub_liao13',tt,ini_x);

figure;
plot(tt,num(:,1),'r.','MarkerSize',15);
% hold on;
% plot(num(:,1),num(:,3),'k.','MarkerSize',15);
% hold on;
% plot(num(:,1),num(:,2),'k.','LineWidth',1);
% hold on;
% plot(num(:,3),num(:,2),'k.','LineWidth',1);
h1=legend('$$x$$');
set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
% tt=0:0.1:200;
% [t,num]=ode45('sub_liao13',tt,ini_x);
% figure;
% plot(tt,num(1:end,1),'r.','MarkerSize',15);
% hold on;
% plot(tt,num(1:end,3),'k.','MarkerSize',15);
% h1=legend('$$h$$','$$\alpha$$','$$\beta$$');
% set(h1,'Interpreter','latex','FontSize',15);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);




