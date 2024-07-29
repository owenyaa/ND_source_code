clear;clc;close all;
load 'xx1.mat';
%% 绘制频谱图动画
for i=1:iii-1
    parameter_a=every_a(i).parameter_a;
    for j=1:N_harm
        amplitude(j,i)=sqrt(parameter_a(j,1)^2+parameter_a(j,2)^2);
    end
    w0=parameter_a(1,1);
    w=1:2:2*N_harm-1;w=w*w0;w=[w;w];w=reshape(w,2*N_harm,1);
    ww(:,i)=w;
    amplitude_w=[zeros(1,N_harm);amplitude(:,i)'];amplitude_w1(:,i)=reshape(amplitude_w,2*N_harm,1);
end
for j=1:iii-1
    for i=1:N_harm
        plot(ww(2*i-1:2*i,j),amplitude_w1(2*i-1:2*i,j),'r.-','LineWidth',1,'MarkerSize',12);
        hold on
    end
    M(j) = getframe;
    cla(1);
end
movie(M,1,3);
writerObj=VideoWriter('test.avi');  %// 定义一个视频文件用来存动画
open(writerObj);
%frame = getframe;            %// 把图像存入视频文件中
writeVideo(writerObj,M);
close(writerObj);
    
%% 绘制相图和时程图动画
Tdata=0:0.01:100;
for k=1:iii-1
    parameter_a=every_a(k).parameter_a;
    x=zeros(N_dof,length(Tdata));dx=zeros(N_dof,length(Tdata));ddx=zeros(N_dof,length(Tdata));
    w0=parameter_a(1,1);
    Harm_parameter_a=parameter_a(2:end,:);
    for j=1:N_dof
        for i=1:N_harm   % i=1,3,5
            x(j,:)=x(j,:)+Harm_parameter_a(i,2*j-1)*cos((2*i-1)*w0*Tdata)+Harm_parameter_a(i,2*j)*sin((2*i-1)*w0*Tdata);
            dx(j,:)=dx(j,:)-w0*(2*i-1)*Harm_parameter_a(i,2*j-1)*sin((2*i-1)*w0*Tdata)+w0*(2*i-1)*Harm_parameter_a(i,2*j)*cos((2*i-1)*w0*Tdata);
            ddx(j,:)=ddx(j,:)-(w0*(2*i-1))^2*Harm_parameter_a(i,2*j-1)*cos((2*i-1)*w0*Tdata)-(w0*(2*i-1))^2*Harm_parameter_a(i,2*j)*sin((2*i-1)*w0*Tdata);
        end
    end
    % 相图
    plot(Tdata,x(1,:),'r-','LineWidth',1);
    hold on;
    plot(Tdata,x(2,:),'k-','LineWidth',1);
%     hold on;
%     plot(x(3,:),dx(3,:),'b-','LineWidth',1);
    h1=legend('$$h$$','$$\alpha$$');
    set(h1,'Interpreter','latex','FontSize',15);
    set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
    
%     %时程图
%     plot(Tdata,x(1,:),'r-','LineWidth',1);
%     hold on;
%     plot(Tdata,x(2,:),'k-','LineWidth',1);
%     hold on;
%     plot(Tdata,x(3,:),'b-','LineWidth',1);
%     h1=legend('$$h$$','$$\alpha$$','$$\beta$$');
%     set(h1,'Interpreter','latex','FontSize',15);
%     set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
    
    M(k) = getframe;
    cla(1);
end
movie(M,1,3);
writerObj=VideoWriter('test2.avi');  %// 定义一个视频文件用来存动画
open(writerObj);
%frame = getframe;            %// 把图像存入视频文件中
writeVideo(writerObj,M);
close(writerObj);

