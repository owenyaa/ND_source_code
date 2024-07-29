clear;clc;close all;
global tf h r1 r4 N_dof N_harm Tdata
global M C K f w parameter_a
N_dof=5;N_harm=10;
r1=2;r4=5;f=1.5;w=0.1;
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
jj=1;
for w=0.1:0.01:10
    tf=2*pi/w;h=2*pi/(w*1000);Tdata=(0:h:tf);
    %% 第一行存储频率，后面存储谐波系数,每个自由度两列
    % parameter_a=[w0,  0,    0,   0;
    %             C_11,S_11,C_21,S_21;
    %             C_12,S_12,C_22,S_22;
    %             C_13,S_13,C_23,S_23;...];
    parameter_a=zeros(N_harm,2*N_dof);
    parameter_a(1,:)=[0,0.03,0,0.01,0,0.06,0,0.05,0,0.07];
    ini_parameter_a=parameter_a;
    iteration=length(Tdata);
    % to indicate where a is in the parametric space
    parameter_a_judge=@(parameter_a)(abs(parameter_a)<10);
    % load simple_fre_data.mat; % load observed data --- Tdata and Xdata
    gammaT=1.414;rhob=0.5; % parameter for trust-region algorithm
    parameter_a_record=parameter_a; % To record the values of parameters during iteration, and for each iteration, a is recorded in a single row of a_record
    TR_record=[];  % recording the parameters during trust region
    %% 本算例采用速度和加速度响应识别，且三个自由度噪声等级不同，线加速度1%，角加速度2%和5%
    %% Response sensitivity iteration
    Nmax=1000;   % maximum number for response sensitivity iteration
    Ntr=20;      % maximum number for trust region iteration
    %% response sensitivity Solution by ode45
    % NT=length(residual(:,1));
    for iii=1:Nmax
        %% 重新载入w0,此处为一个周期内均等选取1K个点
        % compute response and response sensitivity for each incremental
        Etol=1e-10;  % Relative error tolerance for convergence of the RS algorithm
        %%
        residual_iden=cal_residual(parameter_a);
        %% SSS为位移响应灵敏度矩阵，第一列和第二列为残差响应，频率灵敏度从第三列开始
        SSS=reshape(residual_iden(:,N_dof+1:2*N_dof),N_dof*length(Tdata),1);%w0
        for i=1:2*N_harm*N_dof-1
            SSS=[SSS,reshape(residual_iden(:,N_dof*(i+1)+1:N_dof*(i+2)),N_dof*length(Tdata),1)];
        end
        dR=-reshape(residual_iden(:,1:N_dof),N_dof*length(Tdata),1);
        [U,s,V]=csvd(SSS);
        lambda_inverse=l_curve(U,s,dR);
        atemp=parameter_a;
        % trust-region algorithm
        for trust=1:Ntr
            %% 计算的da排序，为w0,C_11,S_11(缺,不参与迭代),C_12,S_12,按照参数矩阵的形式排列da
            real_da=tikhonov(U,s,V,dR,lambda_inverse);
            da=reshape(real_da,2,N_dof*N_harm);da=da';
            sensitivity_parameter_da=da(1:N_harm,1:2);
            for num_dof=1:N_dof-1
                sensitivity_parameter_da=[sensitivity_parameter_da,da(num_dof*N_harm+1:(num_dof+1)*N_harm,1:2)];%da(N_harm+1:2*N_harm,1:2);
            end
            %%
            if ~parameter_a_judge(atemp+sensitivity_parameter_da)      % if updated a is not out of the parametric space, then, lambda should be increased until updated a is in ...
                lambda_inverse=lambda_inverse*gammaT; %  update of lambda
                continue;
            end
            %% 计算的da排序，为w0,C_11,S_11(缺,不参与迭代),C_12,S_12,按照参数矩阵的形式排列da
            real_da=tikhonov(U,s,V,dR,lambda_inverse);
            ini_da=real_da;
            da=reshape(real_da,2,N_dof*N_harm);da=da';
            sensitivity_parameter_da=da(1:N_harm,1:2);
            for num_dof=1:N_dof-1
                sensitivity_parameter_da=[sensitivity_parameter_da,da(num_dof*N_harm+1:(num_dof+1)*N_harm,1:2)];%da(N_harm+1:2*N_harm,1:2);
            end
            %%
            parameter_a=atemp+sensitivity_parameter_da;
            %% 重新计算Tdata
            %% 用新的parameter_a=atemp+da计算响应
            residual_da=cal_residual(parameter_a);
            
            dRtemp=-reshape(residual_da(:,1:N_dof),N_dof*length(Tdata),1);
            LdR=SSS*ini_da-dR;
            rhos=(dR'*dR-dRtemp'*dRtemp)/(dR'*dR-LdR'*LdR);  % agreement indicator
            if rhos>=rhob
                break;
            end
            lambda_inverse=lambda_inverse*gammaT;
        end
        tolt=norm(da)/norm(parameter_a)
        parameter_a_record=[parameter_a_record,parameter_a];
        TR_record=[TR_record;lambda_inverse];
        parameter_a;
        if tolt<=Etol
            break;
        end
        every_a(iii).parameter_a=parameter_a;
        iii
    end
    w
    parameter_a
    bifurcation_p(jj).w=w;
    bifurcation_p(jj).parameter_a=parameter_a;
    jj=jj+1
    residual=cal_residual(parameter_a);
    
    figure;
    plot(Tdata,residual(:,1),'r-','LineWidth',1.5);
    hold on;
    plot(Tdata,residual(:,2),'b-','LineWidth',1.5);
    hold on;
    plot(Tdata,residual(:,3),'k-','LineWidth',1.5);
    hold on;
    plot(Tdata,residual(:,4),'g-','LineWidth',1.5);
    hold on;
    plot(Tdata,residual(:,5),'p-','LineWidth',1.5);
    h1=legend('$$x_1$$','$$x_2$$','$$x_3$$','$$x_4$$','$$x_5$$');
    set(h1,'Interpreter','latex','FontSize',15);
    set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
    
    % Tdata1=0:0.01:150;
    % Harm_parameter_a=parameter_a(1:end,:);
    % %% 计算方程残差
    % x=zeros(N_dof,length(Tdata1));dx=zeros(N_dof,length(Tdata1));ddx=zeros(N_dof,length(Tdata1));
    % for j=1:N_dof
    %     for i=1:N_harm   % i=1,3,5
    %         x(j,:)=x(j,:)+Harm_parameter_a(i,2*j-1)*cos((2*i-1)*w*Tdata1)+Harm_parameter_a(i,2*j)*sin((2*i-1)*w*Tdata1);
    %         dx(j,:)=dx(j,:)-w*(2*i-1)*Harm_parameter_a(i,2*j-1)*sin((2*i-1)*w*Tdata1)+w*(2*i-1)*Harm_parameter_a(i,2*j)*cos((2*i-1)*w*Tdata1);
    %         ddx(j,:)=ddx(j,:)-(w*(2*i-1))^2*Harm_parameter_a(i,2*j-1)*cos((2*i-1)*w*Tdata1)-(w*(2*i-1))^2*Harm_parameter_a(i,2*j)*sin((2*i-1)*w*Tdata1);
    %     end
    % end
    %
    % figure;
    % plot(Tdata1,x(1,:),'r-','LineWidth',1.5);
    % hold on;
    % plot(Tdata1,x(2,:),'b-','LineWidth',1.5);
    % hold on;
    % plot(Tdata1,x(3,:),'k-','LineWidth',1.5);
    % hold on;
    % plot(Tdata1,x(4,:),'g-','LineWidth',1.5);
    % hold on;
    % plot(Tdata1,x(5,:),'p-','LineWidth',1.5);
    % h1=legend('$$x_1$$','$$x_2$$','$$x_3$$','$$x_4$$','$$x_5$$');
    % set(h1,'Interpreter','latex','FontSize',15);
    % set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
    % figure;
    % plot(x(1,:),dx(1,:),'r-','LineWidth',1.5);
    % hold on;
    % plot(x(2,:),dx(2,:),'b-','LineWidth',1.5);
    % hold on;
    % plot(x(3,:),dx(3,:),'k-','LineWidth',1.5);
    % hold on;
    % plot(x(4,:),dx(4,:),'g-','LineWidth',1.5);
    % hold on;
    % plot(x(5,:),dx(5,:),'p-','LineWidth',1.5);
    % h1=legend('$$x_1$$','$$x_2$$','$$x_3$$','$$x_4$$','$$x_5$$');
    % set(h1,'Interpreter','latex','FontSize',15);
    % set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
    
    ini_parameter_a;
end