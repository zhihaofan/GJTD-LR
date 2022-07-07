function [SigmaX,svp]=ClosedWNNM(SigmaY,C,oureps)           % SigmaY:svd特征值矩阵
temp=(SigmaY-oureps).^2-4*(C-oureps*SigmaY);  % .*:对应元素相乘，值为c2
ind=find (temp>0);   % 找到矩阵>0的位置
svp=length(ind);
SigmaX=max(SigmaY(ind)-oureps+sqrt(temp(ind)),0)/2;             % 软阈值

end