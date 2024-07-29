maxCycle=1000; 
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
        %         FitnessSol=fitness1(sol(i,1),sol(i,2));
        [parameter_a]=twoten(sol(i,:));
        FitnessSol=fitness1(parameter_a);
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
                Foods(i,:)=suiji;
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
            end;
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
            %              Fitness(i)=fitness1(sol(t,1),sol(t,2));
            [parameter_a]=twoten(sol(t,:));
            Fitness(i)=fitness1(parameter_a);
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
                    Foods(t,:)=suiji;   %超过设定的Limit次，则该观察蜂变成侦查蜂
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
        Foods(i,:)=suiji;
        %         Fitness(i)=fitness1(Foods(i,1),Foods(i,2)); %求侦查蜂坐标位置的函数值（适应值）
        [parameter_a]=twoten(Foods(i,:));
        Fitness(i)=fitness1(parameter_a);
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
end 
toc   %计时结束
%%%画图%%%
figure(1)
plot(1:maxCycle,Zuiyou)
%axis([0,iter,o,1.2]) %图坐标范围
% title('函数 f(x,y)=(sin(x)/x).*(sin(y)/y)')
xlabel('迭代次数')
ylabel('函数值')
