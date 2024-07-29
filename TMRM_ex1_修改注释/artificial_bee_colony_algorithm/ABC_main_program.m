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
%% 第一行存储频率，后面存储谐波系数,每个自由度两列
parameter_a=zeros(order+1,1);
parameter_a(1,1)=1;parameter_a(2,1)=-0.1;parameter_a(3,1)=0.1;
ini_parameter_a=parameter_a;
%%
% tic                      %计时开始
%%%初始化种群参数%%%
global Optimalvalue;   %定义全局变量，存放寻优过程最优值
Optimalvalue=0;        %给最优变量初值，该变量用以存放所求的最优值
NP=200;                 %初始化蜜蜂总数（采蜜蜂+观察蜂+侦查蜂）
FoodNumber=NP/2;       %确定采蜜蜂数量
SearchNumber=15;        %侦查蜂数量
maxCycle=6000;           %程序总最大迭代次数
limit=10;
CR=0.8;               %DE算法交叉概率
%此处参数为谐波系数的范围，更改generate_random_solution文件是，此处要更改
% a1=1.5;a2=-15;a3=100
% a1_min=0.1;a1_max=3;a2_min=-20;a2_max=-2;a3_min=10;a3_max=150;%ex1
% a1_min=0.1;a1_max=3;a2_min=-50;a2_max=-2;a3_min=10;a3_max=150;%ex2
a_min=-3;a_max=3;
%对每个蜂适应值的不变性进行统计，达到Limit次后进行相关变异操作
%%
%% Author : GUANG_LIU  * owenyaa@gmail.com *
% Created Time : 2016-06-07 13:27
% Last Revised : GUANG_LIU ,2016-06-07
x=order+1;                %自变量数目
Foods=zeros(NP,x);
for i=1:NP
    Foods(i,:)=generate_random_solution;
end
%所有（行）个蜂位置（x,y）坐标矩阵随机值
sol=Foods;                   %定义变量
trial=zeros(NP);             %对各蜂适应值的不变性次数进行统计的变量,zeros产生全为0的方阵
Foodspath=Foods;             %变量保存所有蜂寻优过程，用以画图显示用
%%
for i=1:NP
    %twoto函数作用，将Foods的字符转变为四个参数，其中x为两列矩阵
    [parameter_a]=convert_to_ten(Foods(i,:));
    %fitness函数根据四个参数计算适应值
    Fitness(i)=calculate_fitness(parameter_a);
    %调用函数求各蜂坐标位置的函数值（适应值）
end
%%
Indfit=sort(Fitness);        %对各蜂适应值由小到大排序，sort为从小到大排序
Optimalvalue=Indfit(1);    %保存最优值
iter=1;   %初始化迭代次数
j=1;      %用以初始化结果
Zuiyou=zeros(maxCycle,1);time_record=zeros(maxCycle,1);
toc;
time_inital=toc;
while (iter<=maxCycle)
    tic;
    %%%采蜜蜂模式,在每一个采蜜蜂周围（邻域）寻找更优秀的解
    for i=1:(FoodNumber)
        %采蜜模式，首先选择一列（列数即问题的变量数）
        ParamxChange=fix(rand*x)+1;   %随机选任一列作为变化列，fix函数为向0靠拢取整
        %在采蜜蜂中选择一个蜜蜂
        neighbour=fix(rand*FoodNumber)+1;    %随机选领域蜂
        while (neighbour==i)
            %如果neighbour==i，则再进行一次随机取值，以使两者不相等（不同的蜂）
            neighbour=fix(rand*(FoodNumber))+1;
        end
        sol(i,:)=Foods(i,:);
        %%%在附近区域进行搜索，即改变neighbour行的ParamxChange列中的变量，改变规则为new_x（i，ParamxChange）=x（i，ParamxChange）+rand*(x（neighbour，ParamxChange）-x（i，ParamxChange）)
        sol(i,ParamxChange)=Foods(i,ParamxChange)+(Foods(i,ParamxChange)-Foods(neighbour,ParamxChange))*(rand-0.5)*2;
        %确保新产生的解（蜜源）位置在参数范围内
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
        %计算新蜜源的适应度函数值,倘若领域的更好，则更新采蜜蜂，局部寻优非常棒
        %         FitnessSol=calculate_fitness(sol(i,1),sol(i,2));
        [parameter_a]=convert_to_ten(sol(i,:));
        FitnessSol=calculate_fitness(parameter_a);
        %%
        %使用贪婪准则，寻找最优蜜源
        if (FitnessSol<Fitness(i))
            Foods(i,:)=sol(i,:);
            Fitness(i)=FitnessSol;
            trial(i)=0;
        else
            trial(i)=trial(i)+1;   %若找不到更好的蜜源，则次数加1
            %%
            if trial(i)>limit
                Foods(i,:)=generate_random_solution;
                %超过设定的Limit次，则该采蜜蜂变成为侦查蜂
            end
        end
    end   %End for
    %此次不用轮盘赌法，改用锦标赛机制；锦标赛机制仍旧依赖于轮盘赌法
    probfit=zeros(1,NP);
    for i=1:NP
        for j=1:NP
            if(Fitness(i)>Fitness(j))
                probfit(i)=probfit(i)+1;
            end
        end
    end
    prob=(0.9*probfit/max(probfit));
    %%%根据适应度值计算采蜜蜂被观察蜂跟随的概率%%%
    %轮盘赌法决定观察蜂是否变为采蜜蜂，但是概率设定很好玩
    %此处轮盘机制，每个适应值都有一个概率，最大适应值概率为1，最小至少为0.1，这个比累加法的轮盘要好，此法可以确保每次迭代最优解一定进入下一个循环
    %此处对于最小值问题要修改，最小值问题，则适应度越小越应该被选中
    %     prob=(0.9*Fitness./max(Fitness))+0.1;%此处为最大值问题对应的概率
    %%%prob=(0.9*Fitness./max(Fitness));%%%此处为常规轮盘赌法
    %%%观察蜂%%%
    %处理采蜜蜂以外的蜂，一个采蜜蜂最多招募一个观察蜂，以增加寻优概率，注意，观察蜂数量=总数-采蜜蜂-侦察蜂
    t=FoodNumber;
    i=1;
    while(t<(NP-SearchNumber))%判断条件，只要采蜜蜂小于总蜂-侦察蜂
        t=t+1;%对观察蜂挨个操作
        %注意，此处操作是采蜜蜂才有的操作，即局部寻优，而观察蜂现在有一定几率也可以进行此操作，观察蜂转变的概率就是轮盘赌法
        if(rand>prob(i))   %按概率选择要跟随的采蜜蜂，此处为最小值问题，故大于号；最大值问题小于号
            ParamxChange=fix(rand*x)+1;   %随机选任一列作为变化列
            %neighbour>N/2,确保是对观察蜂进行操作
            neighbour=FoodNumber+fix(rand*(NP-FoodNumber-SearchNumber))+1;
            %选择任意一个观察蜂
            while(neighbour==t)
                %若相等则再观察蜂中进行一次随机取值，这使得neighbour与t不相等
                neighbour=FoodNumber+fix(rand*(NP-FoodNumber-SearchNumber))+1;
            end
            sol(t,:)=Foods(t,:);
            %观察蜂变为对应的采蜜蜂，若变异更好，则留下
            sol(t,ParamxChange)=Foods(i,ParamxChange)+(Foods(i,ParamxChange)-Foods(neighbour,ParamxChange))*(rand-0.5)*2;
            %确保新产生的解在范围内
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
            %计算新蜜源的适应度函数值
            %              Fitness(i)=calculate_fitness(sol(t,1),sol(t,2));
            [parameter_a]=convert_to_ten(sol(t,:));
            Fitness(i)=calculate_fitness(parameter_a);
           %%         
            %使用贪婪准则寻找最优
            if(FitnessSol<Fitness(t))   %若找到更好的蜜源，搜索次数清零
                Foods(t,:)=sol(t,:);
                Fitness(t)=FitnessSol;
                trial(t)=0;
            else
                trial(t)=trial(t)+1;   %若找不到更好的蜜源，则次数加1
                if trial(t)>limit
                    %%
                    Foods(t,:)=generate_random_solution;   %超过设定的Limit次，则该观察蜂变成侦查蜂
                    %%
                end
            end
        end   %End if(rand<prob(i))
        i=i+1;
        if i==(FoodNumber+1)   %当观察蜂数量（t循环次数）远大于采蜜蜂数时，采蜜蜂可能招到多个观察蜂
            i=1;
        end
    end   %End while(t<NP)
    %%%侦查蜂%%%
    %计算侦查蜂位置适应值
    for i=(NP-SearchNumber+1):NP
        %%
        Foods(i,:)=generate_random_solution;
        %         Fitness(i)=calculate_fitness(Foods(i,1),Foods(i,2)); %求侦查蜂坐标位置的函数值（适应值）
        [parameter_a]=convert_to_ten(Foods(i,:));
        Fitness(i)=calculate_fitness(parameter_a);
        %%
    end
    %%%
    %每次循环最后找出最大值
    ind=find(Fitness==min(Fitness));   %找最大值位置
    Min=Fitness(ind);                  %保存最优值
    if (Min<Optimalvalue)
        Optimalvalue=Min;              %保存全局最优值
        GlabolParams=Foods(ind,:);     %保存全局最优值位置
    end
    %%%%%%
    fprintf('iteration=%d,Optimalvalue=%d\n',iter,Optimalvalue);   %输出迭代次数和该次迭代的最优解
    %建立一个迭代次数和最优解的矩阵，以便于画图
    Zuiyou(iter)=Optimalvalue;   %每次迭代得出的最优数值矩阵
    iter=iter+1;
    toc;
    time_record(iter)=toc;
    Foodspath=[Foodspath Foods]; %保存所有蜂寻优位置移动过程
end   %End while
timesum=sum(time_record)+time_inital
% toc   %计时结束
%%%画图%%%
figure;
plot(1:maxCycle,Zuiyou,'r-','LineWidth',1)
%axis([0,iter,o,1.2]) %图坐标范围
% title('函数 f(x,y)=(sin(x)/x).*(sin(y)/y)')
xlabel('迭代次数')
ylabel('函数值')

