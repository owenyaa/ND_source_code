%% Author : GUANG_LIU  * owenyaa@gmail.com *
% Created Time : 2020-11-12 21:08
% Last Revised : GUANG_LIU ,2020-11-12
% Remark : The structure, function and variables of this program are explained:%
%          This document is about the verification procedure of sensitivity matrix,%
%          which verifies the sensitivity response with respect to the parameters by difference.%
clc;clear;close all;
global order Tdata
order=5;Tdata=0:0.01:20;
parameter_a=(0:1:order)';
parameter_a(1,1)=1;
%% ����в�
residual=cal_residual(parameter_a);
%% ���Ʋв�����
% figure;
% plot(Tdata,residual(:,1),'r-','LineWidth',1);
% h1=legend('$$h$$');
% set(h1,'Interpreter','latex','FontSize',15);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);

%% �˴�Ϊ��������֤���򣬱���Ϊ�ò�ִ��������ȣ�������֤�������Ƿ����
% ����������֮�󣬽�ĳ��������ȥһ��С�������²�������һ��
% �������ý�������ٳ���С����Ϊ�����ȣ���֤��ֵ������Ⱥ�ֱ�Ӽ���������������Ƿ��غ�
% ����Ҫע�⣬��ֽ���϶��ǶԵ�(�����Ĳ������)�����ܳ������ԭ��������ݣ�x_cal������
% ������ڶ�Ӧ����ֽ��(parameter_a1)�е�λ��(x1(1,:))���ٶȣ����ٶȶ�Ӧ��ԭ����(parameter_a)�е�����������x_cal(2,:)
% % 
% for i=1:order
%     ddt=0.000001;
%     %��ʼ������
%     parameter_a=(0:1:order)';
%     parameter_a(1,1)=1;
%     parameter_a(i+1,1)=parameter_a(i+1,1)+ddt;% ��������ȥ��
%     
%     residual1=cal_residual(parameter_a);
%     aaaa=(residual1-residual)/ddt;
% %     figure; 
% %     plot(Tdata,residual(:,2*i+2),'r-*')
% %     hold on
% %     plot(Tdata,aaaa(:,2),'k-*')
%     figure; 
%     plot(Tdata,residual(:,2*i+1),'r-*')
%     hold on
%     plot(Tdata,aaaa(:,1),'k-*')
% end


for i=1
    % i=8;
    ddt=0.0001;
    parameter_a=(0:1:order)';
    parameter_a(1,1)=1+ddt;
    residual1=cal_residual(parameter_a);
    aaaa=(residual1-residual)/ddt;
    figure; 
    plot(Tdata,residual(:,end-1),'r-')
    hold on
    plot(Tdata,aaaa(:,1),'k-')
    
end














