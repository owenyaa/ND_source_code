clear;close all;clc;
load one_q10_fu_0_15_fu_0_1.mat;
ii=1;jj=1;
for i=1:length(q10)
    for j=1:length(q20)
        if index(i,j)==1%代表大极限环
            xx(ii)=q10(i);
            yy(ii)=q20(j);
            ii=ii+1;
        else
            xx2(jj)=q10(i);
            yy2(jj)=q20(j);
            jj=jj+1;
        end
    end
end
figure;
plot(xx(1:1:end),yy(1:1:end),'r.');
hold on;
plot(xx2(1:1:end),yy2(1:1:end),'k.');

clear;
load two_q10_fu_0_1_fu_0_05.mat;
ii=1;jj=1;
for i=1:length(q10)
    for j=1:length(q20)
        if index(i,j)==1%代表大极限环
            xx(ii)=q10(i);
            yy(ii)=q20(j);
            ii=ii+1;
        else
            xx2(jj)=q10(i);
            yy2(jj)=q20(j);
            jj=jj+1;
        end
    end
end
hold on;
plot(xx(1:1:end),yy(1:1:end),'r.');
hold on;
plot(xx2(1:1:end),yy2(1:1:end),'k.');

clear;
load three_q10_fu_0_05_0.mat;
ii=1;jj=1;
for i=1:length(q10)
    for j=1:length(q20)
        if index(i,j)==1%代表大极限环
            xx(ii)=q10(i);
            yy(ii)=q20(j);
            ii=ii+1;
        else
            xx2(jj)=q10(i);
            yy2(jj)=q20(j);
            jj=jj+1;
        end
    end
end
hold on;
plot(xx(1:1:end),yy(1:1:end),'r.');
hold on;
plot(xx2(1:1:end),yy2(1:1:end),'k.');

clear;
load four_q10_0_zheng_0_05.mat;
ii=1;jj=1;
for i=1:length(q10)
    for j=1:length(q20)
        if index(i,j)==1%代表大极限环
            xx(ii)=q10(i);
            yy(ii)=q20(j);
            ii=ii+1;
        else
            xx2(jj)=q10(i);
            yy2(jj)=q20(j);
            jj=jj+1;
        end
    end
end
hold on;
plot(xx(1:1:end),yy(1:1:end),'r.');
hold on;
plot(xx2(1:1:end),yy2(1:1:end),'k.');

clear;
load five_q10_zheng_0_05_0_1.mat;
ii=1;jj=1;
for i=1:length(q10)
    for j=1:length(q20)
        if index(i,j)==1%代表大极限环
            xx(ii)=q10(i);
            yy(ii)=q20(j);
            ii=ii+1;
        else
            xx2(jj)=q10(i);
            yy2(jj)=q20(j);
            jj=jj+1;
        end
    end
end
hold on;
plot(xx(1:1:end),yy(1:1:end),'r.');
hold on;
plot(xx2(1:1:end),yy2(1:1:end),'k.');

clear;
load six_q10_zheng_0_1_0_15.mat;
ii=1;jj=1;
for i=1:length(q10)
    for j=1:length(q20)
        if index(i,j)==1%代表大极限环
            xx(ii)=q10(i);
            yy(ii)=q20(j);
            ii=ii+1;
        else
            xx2(jj)=q10(i);
            yy2(jj)=q20(j);
            jj=jj+1;
        end
    end
end
hold on;
plot(xx(1:1:end),yy(1:1:end),'r.');
hold on;
plot(xx2(1:1:end),yy2(1:1:end),'k.');




