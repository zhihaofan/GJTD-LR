
function  [X] =  WNNM( Y, C, NSig)
[U,SigmaY,V] = svd(full(Y),'econ');    
PatNum       = size(Y,2);
%     TempC  = C*sqrt(PatNum)*2*NSig^2;
TempC  = C*NSig^2;
[SigmaX,svp] = ClosedWNNM(SigmaY,TempC,eps); % epsMATLABĬ�ϵ���С����������.    sigmaX��X������ֵ���ֱ�Ϊd1,d2������svp�Ǿ����д���0�ĸ���
X =  U(:,1:svp)*diag(SigmaX)*V(:,1:svp)';     % ������ֵ�ָ�X
return;
