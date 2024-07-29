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
%% 计算残差
residual=cal_residual(parameter_a);
%% 绘制残差曲线
% figure;
% plot(Tdata,residual(:,1),'r-','LineWidth',1);
% h1=legend('$$h$$');
% set(h1,'Interpreter','latex','FontSize',15);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);

%% 此处为灵敏度验证程序，本质为用差分代替灵敏度，可以验证灵敏度是否算错
% 计算正问题之后，将某个参数减去一个小量，用新参数再算一次
% 两次所得结果做差再除以小量即为灵敏度，验证差分的灵敏度和直接计算的灵敏度曲线是否重合
% 但是要注意，差分结果肯定是对的(即下文差量结果)，可能出错的是原程序的内容，x_cal中内容
% 另外关于对应，差分结果(parameter_a1)中的位移(x1(1,:))，速度，加速度对应到原程序(parameter_a)中的灵敏度内容x_cal(2,:)
% % 
% for i=1:order
%     ddt=0.000001;
%     %初始化参数
%     parameter_a=(0:1:order)';
%     parameter_a(1,1)=1;
%     parameter_a(i+1,1)=parameter_a(i+1,1)+ddt;% 增量加上去了
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














