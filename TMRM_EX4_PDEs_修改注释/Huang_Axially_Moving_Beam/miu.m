%% 选择L函数
function l=miu(number,l2,l1,beta)
%g=round(g); % 取四舍五入
% g=ceil(g); % 取g的整数部分+1（若g=4.5,取5）
% g=floor(g); % 取g的整数部分
%number 步数
%l2 一开始的单元数
%l1 最后的个数
%beta 控制收敛快慢的参数
l=[l1+(l2-l1)*exp(-beta*(number-1))];
l=floor(l);
