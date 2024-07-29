function residual=cal_residual(parameter_a)
global N_dof N_harm Tdata 
global M C K k12 k13 k22 k23 w f1 f2
%% 计算方程残差
x=zeros(N_dof,length(Tdata));dx=zeros(N_dof,length(Tdata));ddx=zeros(N_dof,length(Tdata));
for j=1:N_dof
    for i=1:N_harm   % i=1,3,5
        x(j,:)=x(j,:)+parameter_a(i,2*j-1)*cos((2*i-1)*w*Tdata)+parameter_a(i,2*j)*sin((2*i-1)*w*Tdata);
        dx(j,:)=dx(j,:)-w*(2*i-1)*parameter_a(i,2*j-1)*sin((2*i-1)*w*Tdata)+w*(2*i-1)*parameter_a(i,2*j)*cos((2*i-1)*w*Tdata);
        ddx(j,:)=ddx(j,:)-(w*(2*i-1))^2*parameter_a(i,2*j-1)*cos((2*i-1)*w*Tdata)-(w*(2*i-1))^2*parameter_a(i,2*j)*sin((2*i-1)*w*Tdata);
    end
end
residual(1:N_dof,:)=M*ddx+C*dx+K*x+[k12*x(1,:).*x(2,:).^2+k13*x(1,:).^3-f1*cos(w*Tdata);k22*x(2,:).*x(1,:).^2+k23*x(2,:).^3-f2*cos(w*Tdata)];
%% 计算谐波系数的灵敏度
for i=1:2*N_harm*N_dof
    sensitivity_parameter_a1=zeros(2*N_harm*N_dof,1);
    x_a=zeros(N_dof,length(Tdata));dx_a=zeros(N_dof,length(Tdata));ddx_a=zeros(N_dof,length(Tdata));
    sensitivity_parameter_a1(i,1)=1;
    sensitivity_parameter_a1=reshape(sensitivity_parameter_a1,2,N_harm*N_dof);sensitivity_parameter_a1=sensitivity_parameter_a1';
    sensitivity_parameter_a=sensitivity_parameter_a1(1:N_harm,1:2);
    for num_dof=1:N_dof-1
        sensitivity_parameter_a=[sensitivity_parameter_a,sensitivity_parameter_a1(num_dof*N_harm+1:(num_dof+1)*N_harm,1:2)];%da(N_harm+1:2*N_harm,1:2);
    end
    for k=1:N_dof
        for j=1:N_harm
            x_a(k,:)=x_a(k,:)+sensitivity_parameter_a(j,2*k-1)*cos((2*j-1)*w*Tdata)+sensitivity_parameter_a(j,2*k)*sin((2*j-1)*w*Tdata);
            dx_a(k,:)=dx_a(k,:)-w*(2*j-1)*sensitivity_parameter_a(j,2*k-1)*sin((2*j-1)*w*Tdata)+w*(2*j-1)*sensitivity_parameter_a(j,2*k)*cos((2*j-1)*w*Tdata);
            ddx_a(k,:)=ddx_a(k,:)-(w*(2*j-1))^2*sensitivity_parameter_a(j,2*k-1)*cos((2*j-1)*w*Tdata)-(w*(2*j-1))^2*sensitivity_parameter_a(j,2*k)*sin((2*j-1)*w*Tdata);
        end
    end
    residual(N_dof*i+1:N_dof*(i+1),:)=M*ddx_a+C*dx_a+K*x_a+[k12*x(2,:).^2.*x_a(1,:)+2*k12*x(2,:).*x(1,:).*x_a(2,:)+3*k13*x(1,:).^2.*x_a(1,:);...
        k22*x(1,:).^2.*x_a(2,:)+2*k22*x(1,:).*x(2,:).*x_a(1,:)+3*k23*x(2,:).^2.*x_a(2,:)]; 
end

residual=residual';
