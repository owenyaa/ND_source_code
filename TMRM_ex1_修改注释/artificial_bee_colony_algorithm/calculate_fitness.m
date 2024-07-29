%适应度函数
function [f] = calculate_fitness(parameter_a)
global order Tdata x dx
a0=parameter_a(1,1);
%% 自由度数目，谐波数目，参数个数2*N_harm*N_dof
%% 计算方程残差
x=zeros(1,length(Tdata));dx=zeros(1,length(Tdata));
x(1,:)=a0;
for i=1:order
    x(1,:)=x(1,:)+parameter_a(i+1,1).*exp(-i.*Tdata);
    dx(1,:)=dx(1,:)-i*parameter_a(i+1,1).*exp(-i.*Tdata);
end
temp_f(1,:)=dx(1,:)+x(1,:).^2-1;
temp_f(2,:)=sum(parameter_a);
f=norm(temp_f);
end
