%% ѡ��L����
function l=miu(number,l2,l1,beta)
%g=round(g); % ȡ��������
% g=ceil(g); % ȡg����������+1����g=4.5,ȡ5��
% g=floor(g); % ȡg����������
%number ����
%l2 һ��ʼ�ĵ�Ԫ��
%l1 ���ĸ���
%beta �������������Ĳ���
l=[l1+(l2-l1)*exp(-beta*(number-1))];
l=floor(l);
