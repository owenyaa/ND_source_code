maxCycle=1000; 
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
        %         FitnessSol=fitness1(sol(i,1),sol(i,2));
        [parameter_a]=twoten(sol(i,:));
        FitnessSol=fitness1(parameter_a);
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
                Foods(i,:)=suiji;
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
            end;
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
            %              Fitness(i)=fitness1(sol(t,1),sol(t,2));
            [parameter_a]=twoten(sol(t,:));
            Fitness(i)=fitness1(parameter_a);
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
                    Foods(t,:)=suiji;   %�����趨��Limit�Σ���ù۲��������
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
        Foods(i,:)=suiji;
        %         Fitness(i)=fitness1(Foods(i,1),Foods(i,2)); %����������λ�õĺ���ֵ����Ӧֵ��
        [parameter_a]=twoten(Foods(i,:));
        Fitness(i)=fitness1(parameter_a);
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
end 
toc   %��ʱ����
%%%��ͼ%%%
figure(1)
plot(1:maxCycle,Zuiyou)
%axis([0,iter,o,1.2]) %ͼ���귶Χ
% title('���� f(x,y)=(sin(x)/x).*(sin(y)/y)')
xlabel('��������')
ylabel('����ֵ')
