%% Author : GUANG_LIU  * owenyaa@gmail.com *
% Created Time : 2020-11-12 21:32
% Last Revised : GUANG_LIU ,2020-11-12
% Remark : The structure, function and variables of this program are explained:%
%          This program is to calculate the example 1 in the paper "A New Semi-%
%          analytical technique for Nonlinear Systems Based on Response Sensitivity%
%           Analysis", that is, the semi-analytical solution of equation %
%          \begin{equation}%
%          	\label{eq4.2}%
%          	\begin{cases}%
%          		\begin{aligned}%
%          		  & \dot{S}(t)+S^2(t)=1, t\geq 0 \\%
%          		  & S(0)=0                       %
%          		\end{aligned}%
%          	\end{cases}%
%          \end{equation}. The structure of the program is as follows:%
%          The basis function in this program is exponential function \chi_i(t)=e^{-it}%
%          The MATLAB codes about the artificial bee_colony (ABC) algorithm%
%          parameter_a :the coefficients%
%          order: the reserved order%
%          Tdata: the duration%



clear;clc;
%%
tic;
% global a1_min a1_max a2_min a2_max a3_min a3_max 
global a_min a_max
global Tdata order N_dof
N_dof=2;
%% different initial values
order=49;Tdata=0:0.01:200;
%% ��һ�д洢Ƶ�ʣ�����洢г��ϵ��,ÿ�����ɶ�����
parameter_a=zeros(order+1,1);
parameter_a(1,1)=1;parameter_a(2,1)=-0.1;parameter_a(3,1)=0.1;
ini_parameter_a=parameter_a;
%%
% tic                      %��ʱ��ʼ
%%%��ʼ����Ⱥ����%%%
global Optimalvalue;   %����ȫ�ֱ��������Ѱ�Ź�������ֵ
Optimalvalue=0;        %�����ű�����ֵ���ñ������Դ�����������ֵ
NP=200;                 %��ʼ���۷����������۷�+�۲��+���䣩
FoodNumber=NP/2;       %ȷ�����۷�����
SearchNumber=15;        %��������
maxCycle=6000;           %����������������
limit=10;
CR=0.8;               %DE�㷨�������
%�˴�����Ϊг��ϵ���ķ�Χ������generate_random_solution�ļ��ǣ��˴�Ҫ����
% a1=1.5;a2=-15;a3=100
% a1_min=0.1;a1_max=3;a2_min=-20;a2_max=-2;a3_min=10;a3_max=150;%ex1
% a1_min=0.1;a1_max=3;a2_min=-50;a2_max=-2;a3_min=10;a3_max=150;%ex2
a_min=-3;a_max=3;
%��ÿ������Ӧֵ�Ĳ����Խ���ͳ�ƣ��ﵽLimit�κ������ر������
%%
%% Author : GUANG_LIU  * owenyaa@gmail.com *
% Created Time : 2016-06-07 13:27
% Last Revised : GUANG_LIU ,2016-06-07
x=order+1;                %�Ա�����Ŀ
Foods=zeros(NP,x);
for i=1:NP
    Foods(i,:)=generate_random_solution;
end
%���У��У�����λ�ã�x,y������������ֵ
sol=Foods;                   %�������
trial=zeros(NP);             %�Ը�����Ӧֵ�Ĳ����Դ�������ͳ�Ƶı���,zeros����ȫΪ0�ķ���
Foodspath=Foods;             %�����������з�Ѱ�Ź��̣����Ի�ͼ��ʾ��
%%
for i=1:NP
    %twoto�������ã���Foods���ַ�ת��Ϊ�ĸ�����������xΪ���о���
    [parameter_a]=convert_to_ten(Foods(i,:));
    %fitness���������ĸ�����������Ӧֵ
    Fitness(i)=calculate_fitness(parameter_a);
    %���ú������������λ�õĺ���ֵ����Ӧֵ��
end
%%
Indfit=sort(Fitness);        %�Ը�����Ӧֵ��С��������sortΪ��С��������
Optimalvalue=Indfit(1);    %��������ֵ
iter=1;   %��ʼ����������
j=1;      %���Գ�ʼ�����
Zuiyou=zeros(maxCycle,1);time_record=zeros(maxCycle,1);
toc;
time_inital=toc;
while (iter<=maxCycle)
    tic;
    %%%���۷�ģʽ,��ÿһ�����۷���Χ������Ѱ�Ҹ�����Ľ�
    for i=1:(FoodNumber)
        %����ģʽ������ѡ��һ�У�����������ı�������
        ParamxChange=fix(rand*x)+1;   %���ѡ��һ����Ϊ�仯�У�fix����Ϊ��0��£ȡ��
        %�ڲ��۷���ѡ��һ���۷�
        neighbour=fix(rand*FoodNumber)+1;    %���ѡ�����
        while (neighbour==i)
            %���neighbour==i�����ٽ���һ�����ȡֵ����ʹ���߲���ȣ���ͬ�ķ䣩
            neighbour=fix(rand*(FoodNumber))+1;
        end
        sol(i,:)=Foods(i,:);
        %%%�ڸ�������������������ı�neighbour�е�ParamxChange���еı������ı����Ϊnew_x��i��ParamxChange��=x��i��ParamxChange��+rand*(x��neighbour��ParamxChange��-x��i��ParamxChange��)
        sol(i,ParamxChange)=Foods(i,ParamxChange)+(Foods(i,ParamxChange)-Foods(neighbour,ParamxChange))*(rand-0.5)*2;
        %ȷ���²����Ľ⣨��Դ��λ���ڲ�����Χ��
        %%
        %         if (ParamxChange==1)
        %             if (sol(i,ParamxChange)>a1_max||sol(i,ParamxChange)<a1_min)
        %                 sol(i,ParamxChange)=Foods(i,ParamxChange);
        %             end
        %         elseif (ParamxChange==2)
        %             if (sol(i,ParamxChange)>a2_max||sol(i,ParamxChange)<a2_min)
        %                 sol(i,ParamxChange)=Foods(i,ParamxChange);
        %             end
        %         elseif (ParamxChange==3)
        %             if (sol(i,ParamxChange)>a3_max||sol(i,ParamxChange)<a3_min)
        %                 sol(i,ParamxChange)=Foods(i,ParamxChange);
        %             end
        %         end
        
        %%
        %��������Դ����Ӧ�Ⱥ���ֵ,��������ĸ��ã�����²��۷䣬�ֲ�Ѱ�ŷǳ���
        %         FitnessSol=calculate_fitness(sol(i,1),sol(i,2));
        [parameter_a]=convert_to_ten(sol(i,:));
        FitnessSol=calculate_fitness(parameter_a);
        %%
        %ʹ��̰��׼��Ѱ��������Դ
        if (FitnessSol<Fitness(i))
            Foods(i,:)=sol(i,:);
            Fitness(i)=FitnessSol;
            trial(i)=0;
        else
            trial(i)=trial(i)+1;   %���Ҳ������õ���Դ���������1
            %%
            if trial(i)>limit
                Foods(i,:)=generate_random_solution;
                %�����趨��Limit�Σ���ò��۷���Ϊ����
            end
        end
    end   %End for
    %�˴β������̶ķ������ý��������ƣ������������Ծ����������̶ķ�
    probfit=zeros(1,NP);
    for i=1:NP
        for j=1:NP
            if(Fitness(i)>Fitness(j))
                probfit(i)=probfit(i)+1;
            end
        end
    end
    prob=(0.9*probfit/max(probfit));
    %%%������Ӧ��ֵ������۷䱻�۲�����ĸ���%%%
    %���̶ķ������۲���Ƿ��Ϊ���۷䣬���Ǹ����趨�ܺ���
    %�˴����̻��ƣ�ÿ����Ӧֵ����һ�����ʣ������Ӧֵ����Ϊ1����С����Ϊ0.1��������ۼӷ�������Ҫ�ã��˷�����ȷ��ÿ�ε������Ž�һ��������һ��ѭ��
    %�˴�������Сֵ����Ҫ�޸ģ���Сֵ���⣬����Ӧ��ԽСԽӦ�ñ�ѡ��
    %     prob=(0.9*Fitness./max(Fitness))+0.1;%�˴�Ϊ���ֵ�����Ӧ�ĸ���
    %%%prob=(0.9*Fitness./max(Fitness));%%%�˴�Ϊ�������̶ķ�
    %%%�۲��%%%
    %������۷�����ķ䣬һ�����۷������ļһ���۲�䣬������Ѱ�Ÿ��ʣ�ע�⣬�۲������=����-���۷�-����
    t=FoodNumber;
    i=1;
    while(t<(NP-SearchNumber))%�ж�������ֻҪ���۷�С���ܷ�-����
        t=t+1;%�Թ۲�䰤������
        %ע�⣬�˴������ǲ��۷���еĲ��������ֲ�Ѱ�ţ����۲��������һ������Ҳ���Խ��д˲������۲��ת��ĸ��ʾ������̶ķ�
        if(rand>prob(i))   %������ѡ��Ҫ����Ĳ��۷䣬�˴�Ϊ��Сֵ���⣬�ʴ��ںţ����ֵ����С�ں�
            ParamxChange=fix(rand*x)+1;   %���ѡ��һ����Ϊ�仯��
            %neighbour>N/2,ȷ���ǶԹ۲����в���
            neighbour=FoodNumber+fix(rand*(NP-FoodNumber-SearchNumber))+1;
            %ѡ������һ���۲��
            while(neighbour==t)
                %��������ٹ۲���н���һ�����ȡֵ����ʹ��neighbour��t�����
                neighbour=FoodNumber+fix(rand*(NP-FoodNumber-SearchNumber))+1;
            end
            sol(t,:)=Foods(t,:);
            %�۲���Ϊ��Ӧ�Ĳ��۷䣬��������ã�������
            sol(t,ParamxChange)=Foods(i,ParamxChange)+(Foods(i,ParamxChange)-Foods(neighbour,ParamxChange))*(rand-0.5)*2;
            %ȷ���²����Ľ��ڷ�Χ��
            %%
            %         if (ParamxChange==1)
            %             if (sol(t,ParamxChange)>a1_max||sol(t,ParamxChange)<a1_min)
            %                 sol(t,ParamxChange)=Foods(t,ParamxChange);
            %             end
            %         elseif (ParamxChange==2)
            %             if (sol(t,ParamxChange)>a2_max||sol(t,ParamxChange)<a2_min)
            %                 sol(t,ParamxChange)=Foods(t,ParamxChange);
            %             end
            %         elseif (ParamxChange==3)
            %             if (sol(t,ParamxChange)>a3_max||sol(t,ParamxChange)<a3_min)
            %                 sol(t,ParamxChange)=Foods(t,ParamxChange);
            %             end
            %         end
            %%
            %��������Դ����Ӧ�Ⱥ���ֵ
            %              Fitness(i)=calculate_fitness(sol(t,1),sol(t,2));
            [parameter_a]=convert_to_ten(sol(t,:));
            Fitness(i)=calculate_fitness(parameter_a);
           %%         
            %ʹ��̰��׼��Ѱ������
            if(FitnessSol<Fitness(t))   %���ҵ����õ���Դ��������������
                Foods(t,:)=sol(t,:);
                Fitness(t)=FitnessSol;
                trial(t)=0;
            else
                trial(t)=trial(t)+1;   %���Ҳ������õ���Դ���������1
                if trial(t)>limit
                    %%
                    Foods(t,:)=generate_random_solution;   %�����趨��Limit�Σ���ù۲��������
                    %%
                end
            end
        end   %End if(rand<prob(i))
        i=i+1;
        if i==(FoodNumber+1)   %���۲��������tѭ��������Զ���ڲ��۷���ʱ�����۷�����е�����۲��
            i=1;
        end
    end   %End while(t<NP)
    %%%����%%%
    %��������λ����Ӧֵ
    for i=(NP-SearchNumber+1):NP
        %%
        Foods(i,:)=generate_random_solution;
        %         Fitness(i)=calculate_fitness(Foods(i,1),Foods(i,2)); %����������λ�õĺ���ֵ����Ӧֵ��
        [parameter_a]=convert_to_ten(Foods(i,:));
        Fitness(i)=calculate_fitness(parameter_a);
        %%
    end
    %%%
    %ÿ��ѭ������ҳ����ֵ
    ind=find(Fitness==min(Fitness));   %�����ֵλ��
    Min=Fitness(ind);                  %��������ֵ
    if (Min<Optimalvalue)
        Optimalvalue=Min;              %����ȫ������ֵ
        GlabolParams=Foods(ind,:);     %����ȫ������ֵλ��
    end
    %%%%%%
    fprintf('iteration=%d,Optimalvalue=%d\n',iter,Optimalvalue);   %������������͸ôε��������Ž�
    %����һ���������������Ž�ľ����Ա��ڻ�ͼ
    Zuiyou(iter)=Optimalvalue;   %ÿ�ε����ó���������ֵ����
    iter=iter+1;
    toc;
    time_record(iter)=toc;
    Foodspath=[Foodspath Foods]; %�������з�Ѱ��λ���ƶ�����
end   %End while
timesum=sum(time_record)+time_inital
% toc   %��ʱ����
%%%��ͼ%%%
figure;
plot(1:maxCycle,Zuiyou,'r-','LineWidth',1)
%axis([0,iter,o,1.2]) %ͼ���귶Χ
% title('���� f(x,y)=(sin(x)/x).*(sin(y)/y)')
xlabel('��������')
ylabel('����ֵ')

