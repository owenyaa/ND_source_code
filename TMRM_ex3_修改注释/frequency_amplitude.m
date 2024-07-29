clear;clc;close all;
load 'N12.mat';
%% »æÖÆÆµÆ×Í¼¶¯»­
for j=1:N_harm
    amplitude1(1,j)=sqrt(parameter_a(j+1,1)^2+parameter_a(j+1,2)^2);
    amplitude2(1,j)=sqrt(parameter_a(j+1,3)^2+parameter_a(j+1,4)^2);
end
w0=parameter_a(1,1);
w=1:2:2*N_harm-1;w=w*w0;w=[w;w];
amplitude_w1(1,:)=zeros(length(amplitude1),1);amplitude_w1(2,:)=log10(amplitude1);amplitude_w1(2,:)=amplitude1;
amplitude_w2(1,:)=zeros(length(amplitude2),1);amplitude_w2(2,:)=log10(amplitude2);amplitude_w2(2,:)=amplitude2;
% amplitude_w1=log10(amplitude_w1);
subplot(2,1,1);
for i=1:N_harm
    %     plot(w,amplitude_w1','r-','LineWidth',1,'MarkerSize',12);
    plot(w(:,i),amplitude_w1(:,i),'r-','LineWidth',1.5);
    hold on
    plot(w(1,i),amplitude_w1(2,i),'k.','MarkerSize',15);
    hold on
end
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
subplot(2,1,2);
for i=1:N_harm
    %     plot(w,amplitude_w1','r-','LineWidth',1,'MarkerSize',12);
    plot(w(:,i),amplitude_w2(:,i),'r-','LineWidth',1.5);
    hold on
    plot(w(1,i),amplitude_w2(2,i),'k.','MarkerSize',15);
    hold on
end

set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);

