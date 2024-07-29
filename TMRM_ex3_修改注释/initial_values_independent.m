clear;clc;close all;
load N12_3.mat;
global tf h Tdata parameter_a beta r N_dof N_harm
%% different initial values
% case 1
N_dof=2;N_harm=12;
beta=20;r=0;
temp_w0=0.01:0.001:1;
U=8;
M=[1,0.25;0.25,0.5];
C=[0.1,0;0,0.1];
K=[0.2,0.1*U;0,0.5-0.04*U];
% for ii=1:length(temp_w0)
%     w0=temp_w0(ii);
%     tf=2*pi/w0;h=2*pi/(w0*1000);Tdata=(0:h:tf);
%     Harm_parameter_a=parameter_a(2:end,:);
%     %% 计算方程残差
%     x=zeros(N_dof,length(Tdata));dx=zeros(N_dof,length(Tdata));ddx=zeros(N_dof,length(Tdata));
%     % 有X,Y两个自由度
%     % for i=1:N_harm   % i=1,3,5
%     for j=1:N_dof
%         for i=1:N_harm   % i=1,3,5
%             x(j,:)=x(j,:)+Harm_parameter_a(i,2*j-1)*cos((2*i-1)*w0*Tdata)+Harm_parameter_a(i,2*j)*sin((2*i-1)*w0*Tdata);
%             dx(j,:)=dx(j,:)-w0*(2*i-1)*Harm_parameter_a(i,2*j-1)*sin((2*i-1)*w0*Tdata)+w0*(2*i-1)*Harm_parameter_a(i,2*j)*cos((2*i-1)*w0*Tdata);
%             ddx(j,:)=ddx(j,:)-(w0*(2*i-1))^2*Harm_parameter_a(i,2*j-1)*cos((2*i-1)*w0*Tdata)-(w0*(2*i-1))^2*Harm_parameter_a(i,2*j)*sin((2*i-1)*w0*Tdata);
%         end
%     end
%     residual(1:N_dof,:)=M*ddx+C*dx+K*x+[r,0;0,beta]*x.^3;
%     %residual=cal_residual(parameter_a);
%     R(ii)=norm(residual);
% end
% 
% 
% temp_a1=-0.5:0.001:1;
% for ii=1:length(temp_a1)
%     w0=parameter_a(1,1);
%     tf=2*pi/w0;h=2*pi/(w0*1000);Tdata=(0:h:tf);
%     Harm_parameter_a=parameter_a(2:end,:);
%     Harm_parameter_a(1,1)=temp_a1(ii);
%     %% 计算方程残差
%     x=zeros(N_dof,length(Tdata));dx=zeros(N_dof,length(Tdata));ddx=zeros(N_dof,length(Tdata));
%     % 有X,Y两个自由度
%     % for i=1:N_harm   % i=1,3,5
%     for j=1:N_dof
%         for i=1:N_harm   % i=1,3,5
%             x(j,:)=x(j,:)+Harm_parameter_a(i,2*j-1)*cos((2*i-1)*w0*Tdata)+Harm_parameter_a(i,2*j)*sin((2*i-1)*w0*Tdata);
%             dx(j,:)=dx(j,:)-w0*(2*i-1)*Harm_parameter_a(i,2*j-1)*sin((2*i-1)*w0*Tdata)+w0*(2*i-1)*Harm_parameter_a(i,2*j)*cos((2*i-1)*w0*Tdata);
%             ddx(j,:)=ddx(j,:)-(w0*(2*i-1))^2*Harm_parameter_a(i,2*j-1)*cos((2*i-1)*w0*Tdata)-(w0*(2*i-1))^2*Harm_parameter_a(i,2*j)*sin((2*i-1)*w0*Tdata);
%         end
%     end
%     residual(1:N_dof,:)=M*ddx+C*dx+K*x+[r,0;0,beta]*x.^3;
%     %residual=cal_residual(parameter_a);
%     R1(ii)=norm(residual);
% end
% 
% temp_b1=-0.5:0.001:1;
% for ii=1:length(temp_b1)
%     w0=parameter_a(1,1);
%     tf=2*pi/w0;h=2*pi/(w0*1000);Tdata=(0:h:tf);
%     Harm_parameter_a=parameter_a(2:end,:);
%     Harm_parameter_a(1,3)=temp_b1(ii);
%     %% 计算方程残差
%     x=zeros(N_dof,length(Tdata));dx=zeros(N_dof,length(Tdata));ddx=zeros(N_dof,length(Tdata));
%     % 有X,Y两个自由度
%     % for i=1:N_harm   % i=1,3,5
%     for j=1:N_dof
%         for i=1:N_harm   % i=1,3,5
%             x(j,:)=x(j,:)+Harm_parameter_a(i,2*j-1)*cos((2*i-1)*w0*Tdata)+Harm_parameter_a(i,2*j)*sin((2*i-1)*w0*Tdata);
%             dx(j,:)=dx(j,:)-w0*(2*i-1)*Harm_parameter_a(i,2*j-1)*sin((2*i-1)*w0*Tdata)+w0*(2*i-1)*Harm_parameter_a(i,2*j)*cos((2*i-1)*w0*Tdata);
%             ddx(j,:)=ddx(j,:)-(w0*(2*i-1))^2*Harm_parameter_a(i,2*j-1)*cos((2*i-1)*w0*Tdata)-(w0*(2*i-1))^2*Harm_parameter_a(i,2*j)*sin((2*i-1)*w0*Tdata);
%         end
%     end
%     residual(1:N_dof,:)=M*ddx+C*dx+K*x+[r,0;0,beta]*x.^3;
%     %residual=cal_residual(parameter_a);
%     R2(ii)=norm(residual);
% end
% 
% temp_c1=-0.5:0.001:1;
% for ii=1:length(temp_c1)
%     w0=parameter_a(1,1);
%     tf=2*pi/w0;h=2*pi/(w0*1000);Tdata=(0:h:tf);
%     Harm_parameter_a=parameter_a(2:end,:);
%     Harm_parameter_a(1,4)=temp_c1(ii);
%     %% 计算方程残差
%     x=zeros(N_dof,length(Tdata));dx=zeros(N_dof,length(Tdata));ddx=zeros(N_dof,length(Tdata));
%     % 有X,Y两个自由度
%     % for i=1:N_harm   % i=1,3,5
%     for j=1:N_dof
%         for i=1:N_harm   % i=1,3,5
%             x(j,:)=x(j,:)+Harm_parameter_a(i,2*j-1)*cos((2*i-1)*w0*Tdata)+Harm_parameter_a(i,2*j)*sin((2*i-1)*w0*Tdata);
%             dx(j,:)=dx(j,:)-w0*(2*i-1)*Harm_parameter_a(i,2*j-1)*sin((2*i-1)*w0*Tdata)+w0*(2*i-1)*Harm_parameter_a(i,2*j)*cos((2*i-1)*w0*Tdata);
%             ddx(j,:)=ddx(j,:)-(w0*(2*i-1))^2*Harm_parameter_a(i,2*j-1)*cos((2*i-1)*w0*Tdata)-(w0*(2*i-1))^2*Harm_parameter_a(i,2*j)*sin((2*i-1)*w0*Tdata);
%         end
%     end
%     residual(1:N_dof,:)=M*ddx+C*dx+K*x+[r,0;0,beta]*x.^3;
%     %residual=cal_residual(parameter_a);
%     R3(ii)=norm(residual);
% end
% [min_w_y,min_w_x]=min(R);
% [min_a1_y,min_a1_x]=min(R1);
% [min_b1_y,min_b1_x]=min(R2);
% [min_c1_y,min_c1_x]=min(R3);
% plot(temp_w0,R,'k-','LineWidth',1.5);
% hold on
% plot(temp_a1,R1,'r-','LineWidth',1.5);
% hold on
% plot(temp_b1,R2,'g-','LineWidth',1.5);
% hold on
% plot(temp_c1,R3,'b-','LineWidth',1.5);
% 
% plot(temp_w0(min_w_x),min_w_y,'k.','MarkerSize',15);
% hold on
% plot(temp_a1(min_a1_x),min_a1_y,'k.','MarkerSize',15);
% hold on
% plot(temp_b1(min_b1_x),min_b1_y,'k.','MarkerSize',15);
% hold on
% plot(temp_c1(min_c1_x),min_c1_y,'k.','MarkerSize',15);
% hold on
% line([-0.5 temp_w0(min_w_x)],[min_w_y min_w_y],'linestyle','-');
% line([temp_w0(min_w_x) temp_w0(min_w_x)],[-0.5 min_w_y],'linestyle','-');
% 
% line([-0.5 temp_a1(min_a1_x)],[min_a1_y min_a1_y],'linestyle','-');
% line([temp_a1(min_a1_x) temp_a1(min_a1_x)],[-0.5 min_a1_y],'linestyle','-');
% 
% line([-0.5 temp_b1(min_b1_x)],[min_b1_y min_b1_y],'linestyle','-');
% line([temp_b1(min_b1_x) temp_b1(min_b1_x)],[-0.5 min_b1_y],'linestyle','-');
% 
% line([-0.5 temp_c1(min_c1_x)],[min_c1_y min_c1_y],'linestyle','-');
% line([temp_c1(min_c1_x) temp_c1(min_c1_x)],[-0.5 min_c1_y],'linestyle','-');
% 
% h1=legend('$$\omega$$','$$b_1$$','$$d_1$$','$$e_1$$');
% set(h1,'Interpreter','latex','FontSize',15);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);

temp_a10=-0.5:0.01:1;
temp_w00=0.01:0.01:1;
for ii=1:length(temp_w00)
    for jj=1:length(temp_a10)
        w0=temp_w00(ii);
        tf=2*pi/w0;h=2*pi/(w0*1000);Tdata=(0:h:tf);
        Harm_parameter_a=parameter_a(2:end,:);
        Harm_parameter_a(1,1)=temp_a10(jj);
        %% 计算方程残差
        x=zeros(N_dof,length(Tdata));dx=zeros(N_dof,length(Tdata));ddx=zeros(N_dof,length(Tdata));
        % 有X,Y两个自由度
        % for i=1:N_harm   % i=1,3,5
        for j=1:N_dof
            for i=1:N_harm   % i=1,3,5
                x(j,:)=x(j,:)+Harm_parameter_a(i,2*j-1)*cos((2*i-1)*w0*Tdata)+Harm_parameter_a(i,2*j)*sin((2*i-1)*w0*Tdata);
                dx(j,:)=dx(j,:)-w0*(2*i-1)*Harm_parameter_a(i,2*j-1)*sin((2*i-1)*w0*Tdata)+w0*(2*i-1)*Harm_parameter_a(i,2*j)*cos((2*i-1)*w0*Tdata);
                ddx(j,:)=ddx(j,:)-(w0*(2*i-1))^2*Harm_parameter_a(i,2*j-1)*cos((2*i-1)*w0*Tdata)-(w0*(2*i-1))^2*Harm_parameter_a(i,2*j)*sin((2*i-1)*w0*Tdata);
            end
        end
        residual(1:N_dof,:)=M*ddx+C*dx+K*x+[r,0;0,beta]*x.^3;
        %residual=cal_residual(parameter_a);
        R10(ii,jj)=norm(residual);
    end
end
figure;
surf(temp_a10,temp_w00,R10)%,'r-','LineWidth',1);
h1=legend('$$a_0$$');
set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);

% [C1,h] = contour(temp_a10,temp_w00,R10);
% set(h,'ShowText','on','TextStep',get(h,'LevelStep')*0.2)
% colormap cool
% 
% for i=1:iii-1
%     w0_record(i)=parameter_a_record(1,4*i-3);
%     a1_record(i)=parameter_a_record(2,4*i-3);
% end
% % hold on
% % plot(a1_record,w0_record,'k-.','MarkerSize',15);
% 
% for ii=1:length(w0_record)
%     for jj=1:length(a1_record)
%         w0=w0_record(ii);
%         tf=2*pi/w0;h=2*pi/(w0*1000);Tdata=(0:h:tf);
%         Harm_parameter_a=parameter_a(2:end,:);
%         Harm_parameter_a(1,1)=a1_record(jj);
%         %% 计算方程残差
%         x=zeros(N_dof,length(Tdata));dx=zeros(N_dof,length(Tdata));ddx=zeros(N_dof,length(Tdata));
%         % 有X,Y两个自由度
%         % for i=1:N_harm   % i=1,3,5
%         for j=1:N_dof
%             for i=1:N_harm   % i=1,3,5
%                 x(j,:)=x(j,:)+Harm_parameter_a(i,2*j-1)*cos((2*i-1)*w0*Tdata)+Harm_parameter_a(i,2*j)*sin((2*i-1)*w0*Tdata);
%                 dx(j,:)=dx(j,:)-w0*(2*i-1)*Harm_parameter_a(i,2*j-1)*sin((2*i-1)*w0*Tdata)+w0*(2*i-1)*Harm_parameter_a(i,2*j)*cos((2*i-1)*w0*Tdata);
%                 ddx(j,:)=ddx(j,:)-(w0*(2*i-1))^2*Harm_parameter_a(i,2*j-1)*cos((2*i-1)*w0*Tdata)-(w0*(2*i-1))^2*Harm_parameter_a(i,2*j)*sin((2*i-1)*w0*Tdata);
%             end
%         end
%         residual(1:N_dof,:)=M*ddx+C*dx+K*x+[r,0;0,beta]*x.^3;
%         %residual=cal_residual(parameter_a);
%         R101(ii,jj)=norm(residual);
%     end
% end
% hold on
% plot3(a1_record(length(a1_record)),w0_record(length(a1_record)),R101(length(a1_record)),'k.','MarkerSize',15);





