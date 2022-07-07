
function  [X] =  WNNM( Y, C, NSig)
[U,SigmaY,V] = svd(full(Y),'econ');    
PatNum       = size(Y,2);
%     TempC  = C*sqrt(PatNum)*2*NSig^2;
TempC  = C*NSig^2;
[SigmaX,svp] = ClosedWNNM(SigmaY,TempC,eps); % epsMATLAB默认的最小浮点数精度.    sigmaX是X的奇异值，分别为d1,d2……，svp是矩阵中大于0的个数
X =  U(:,1:svp)*diag(SigmaX)*V(:,1:svp)';     % 用奇异值恢复X
return;
