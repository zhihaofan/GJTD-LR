function [SigmaX,svp]=ClosedWNNM(SigmaY,C,oureps)           % SigmaY:svd����ֵ����
temp=(SigmaY-oureps).^2-4*(C-oureps*SigmaY);  % .*:��ӦԪ����ˣ�ֵΪc2
ind=find (temp>0);   % �ҵ�����>0��λ��
svp=length(ind);
SigmaX=max(SigmaY(ind)-oureps+sqrt(temp(ind)),0)/2;             % ����ֵ

end