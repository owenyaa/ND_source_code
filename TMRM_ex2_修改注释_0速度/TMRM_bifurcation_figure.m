clear;%close all;
load 'bifurcation_4.2_12.mat';
N_dof=2;N_harm=15;
for ii=1:4:157
    parameter_a=every_a(ii).parameter_a;
    Q(ii)=every_a(ii).Q;
    Tdata=0:0.01:100;
    w0=parameter_a(1,1);
    Harm_parameter_a=parameter_a(2:end,:);
    %% 计算方程残差
    x=zeros(N_dof,length(Tdata));dx=zeros(N_dof,length(Tdata));ddx=zeros(N_dof,length(Tdata));
    % 有X,Y两个自由度
    % for i=1:N_harm   % i=1,3,5
    for j=1:N_dof
        for i=1:N_harm   % i=1,3,5
            x(j,:)=x(j,:)+Harm_parameter_a(i,2*j-1)*cos((2*i-1)*w0*Tdata)+Harm_parameter_a(i,2*j)*sin((2*i-1)*w0*Tdata);
%             dx(j,:)=dx(j,:)-w0*(2*i-1)*Harm_parameter_a(i,2*j-1)*sin((2*i-1)*w0*Tdata)+w0*(2*i-1)*Harm_parameter_a(i,2*j)*cos((2*i-1)*w0*Tdata);
%             ddx(j,:)=ddx(j,:)-(w0*(2*i-1))^2*Harm_parameter_a(i,2*j-1)*cos((2*i-1)*w0*Tdata)-(w0*(2*i-1))^2*Harm_parameter_a(i,2*j)*sin((2*i-1)*w0*Tdata);
        end
    end
    x=x';
    xmax_2=x(:,2);
    xmax(ii).xmax=getmax(xmax_2);
    xmin(ii).xmin=getmin(xmax_2);
    hold on;
    QQ=Q(ii)*ones(1,length(xmax(ii).xmax));
    plot(QQ,xmax(ii).xmax,'r.','MarkerSize',12);
    hold on;
    QQ1=Q(ii)*ones(1,length(xmin(ii).xmin));
    plot(QQ1,xmin(ii).xmin,'r.','MarkerSize',12);
    hold on;
%     figure;
%     plot(Tdata,x(1,:),'k-','LineWidth',1);
%     hold on;
%     plot(Tdata,x(2,:),'b-','LineWidth',1);
    
end
% h1=legend('$$h$$','$$\alpha$$');
% set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);









