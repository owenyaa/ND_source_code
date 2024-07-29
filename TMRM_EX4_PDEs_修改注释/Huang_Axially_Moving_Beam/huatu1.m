clear;
global Tdata parameter_a N_dof N_harm N_w0 index_global Q N_min alpha
%% different initial values
% 不同系统系需要更改的参数
N_dof=3;N_harm=128;N_w0=2;%基频个数
N_min=60; alpha=0.03;
Q=8;
%% 计算基频的组合系数
index_global=[1 0;0 1];
for N=3:2:15
    temp=[N 0];
    for n=N-1:-1:1
        m=N-n;
        temp=[temp;n m;n -m];
    end
    temp=[temp;0 N];
    index_global=[index_global;temp];
end
index_global=[index_global,index_global,index_global];

% load 'no_ini_parameter_a.mat';
load '2_10_7_128.mat';
Tdata=0:0.01:200;
residual=cal_residual(parameter_a);
x=zeros(N_dof,length(Tdata));dx=zeros(N_dof,length(Tdata));ddx=zeros(N_dof,length(Tdata));
% 有X,Y两个自由度
fundamental_w=parameter_a(1,1:N_w0);
vector_w=index_global(:,1:N_w0)*fundamental_w';
Harm_parameter_a=parameter_a(2:end,:);
for j=1:N_dof
    for i=1:N_harm   % i=1,3,5
        x(j,:)=x(j,:)+Harm_parameter_a(i,2*j-1)*cos(vector_w(i)*Tdata)+Harm_parameter_a(i,2*j)*sin(vector_w(i)*Tdata);
        dx(j,:)=dx(j,:)-vector_w(i)*Harm_parameter_a(i,2*j-1)*sin(vector_w(i)*Tdata)+vector_w(i)*Harm_parameter_a(i,2*j)*cos(vector_w(i)*Tdata);
        ddx(j,:)=ddx(j,:)-(vector_w(i))^2*Harm_parameter_a(i,2*j-1)*cos(vector_w(i)*Tdata)-(vector_w(i))^2*Harm_parameter_a(i,2*j)*sin(vector_w(i)*Tdata);
    end
end

for j=1:N_dof
    for i=1:N_harm
        vector_amplitude(i,j)=sqrt(Harm_parameter_a(i,2*j-1)^2+Harm_parameter_a(i,2*j)^2);
    end
end
figure;
plot(abs(vector_w)./(2*pi),vector_amplitude(:,3),'k.','LineWidth',1.5);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);

figure;
plot(Tdata,residual(:,1),'r-','LineWidth',1);
hold on;
plot(Tdata,residual(:,2),'k-','LineWidth',1);
hold on;
plot(Tdata,residual(:,3),'b-','LineWidth',1);
h1=legend('$$h$$','$$\alpha$$','$$\beta$$');
set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);

figure;
subplot(3,1,1);
plot(Tdata,residual(:,1),'k-','LineWidth',1);
subplot(3,1,2);
plot(Tdata,residual(:,2),'k-','LineWidth',1);
subplot(3,1,3);
plot(Tdata,residual(:,3),'k-','LineWidth',1);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);







