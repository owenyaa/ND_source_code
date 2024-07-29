


N_harm=5;
global n_harm%

n_harm=2*N_harm;
% wn=2*pi/period;
period=2;
wn=2*pi/period;
temp1=fix(length(tt)*2/3);temp2=length(tt);
w_input=wn;
t_input=tt(temp1:temp2);
vector=[];
for nx=1:2
    temp=converse(time_to_fourier(w_input,t_input,num(temp1:temp2,nx)));
    vector=[vector temp];%....x_input, for Jacobi and residual
end


% % vector
% vector(1,:)=zeros(size(vector(1,:)));
% temp_el(1)=1;temp_el_i=2;
% for i=1:n_harm
%     if mod(i,2)==0
%         vector(2*i,:)=zeros(size(vector(2*i,:)));
%         vector(2*i+1,:)=zeros(size(vector(2*i+1,:)));
%         temp_el(temp_el_i)=2*i;temp_el(temp_el_i+1)=2*i+1;
%         temp_el_i=temp_el_i+2;
%     end
% end
% % vector
% 
% %开始将vector转换为parameter_a
% parameter_a(1,1)=wn;
% for i=1:2
%     jj=2;
%     for j=1:n_harm
%         if mod(j,2)==1
%             parameter_a(jj,2*i-1)=vector(2*j,i);
%             parameter_a(jj,2*i)=vector(2*j+1,i);
%             jj=jj+1;
%         end
%     end
% end

%开始将vector转换为parameter_a
parameter_a(1,1)=wn;
for i=1:2
    jj=2;
    for j=1:n_harm
        %if mod(j,2)==1
        parameter_a(jj,2*i-1)=vector(2*j,i);
        parameter_a(jj,2*i)=vector(2*j+1,i);
        jj=jj+1;
        %end
    end
end



% N_dof=2;
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


N_dof=2;
Tdata1=0:0.01:150;
Harm_parameter_a=parameter_a(1:end,:);
%% 计算方程残差
x=zeros(N_dof,length(Tdata1));dx=zeros(N_dof,length(Tdata1));ddx=zeros(N_dof,length(Tdata1));
for j=1:N_dof
    for i=1:n_harm   % i=1,3,5
        x(j,:)=x(j,:)+Harm_parameter_a(i,2*j-1)*cos(i*w*Tdata1)+Harm_parameter_a(i,2*j)*sin(i*w*Tdata1);
        dx(j,:)=dx(j,:)-w*i*Harm_parameter_a(i,2*j-1)*sin(i*w*Tdata1)+w*i*Harm_parameter_a(i,2*j)*cos(i*w*Tdata1);
        ddx(j,:)=ddx(j,:)-(w*i)^2*Harm_parameter_a(i,2*j-1)*cos(i*w*Tdata1)-(w*i)^2*Harm_parameter_a(i,2*j)*sin(i*w*Tdata1);
    end
end


figure;
plot(Tdata1,x(1,:),'r-','LineWidth',1.5);
hold on;
plot(Tdata1,x(2,:),'b-','LineWidth',1.5);
h1=legend('$$x_1$$','$$x_2$$','$$x_3$$','$$x_4$$','$$x_5$$');
set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);







