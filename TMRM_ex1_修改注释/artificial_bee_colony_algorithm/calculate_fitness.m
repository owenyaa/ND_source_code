%��Ӧ�Ⱥ���
function [f] = calculate_fitness(parameter_a)
global order Tdata x dx
a0=parameter_a(1,1);
%% ���ɶ���Ŀ��г����Ŀ����������2*N_harm*N_dof
%% ���㷽�̲в�
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
