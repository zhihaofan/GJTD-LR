function [Dict,Drls,coef] = FDDL(TrainDat,TrainLabel,opts)
%%   进度：over
%%%%%%%%%%%%%%%%%%
% normalize energy
%%%%%%%%%%%%%%%%%%
[m1,n1,k1]=size(TrainDat);
for i=1:k1
    TrainDat(:,:,i) = TrainDat(:,:,i)*diag(1./sqrt(sum(TrainDat(:,:,i).*TrainDat(:,:,i))));
end

% TrainDat = TrainDat*diag(1./sqrt(sum(TrainDat.*TrainDat)));

%%%%%%%%%%%%%%%%%%
%initialize dict  
%%%%%%%%%%%%%%%%%% 
Dict_ini  =  [];   %初始化的字典为空，用来存放初始化的字典
Dlabel_ini = [];   %对应的标签

for ci = 1:opts.nClass
    cdat          =    TrainDat(:,TrainLabel==ci,:);    %第一类：300*7 即：7个300维的样本
%     dict          =    FDDL_INID(cdat,size(cdat,2),opts.wayInit);  %对第一类（7个样本）而言，用pca初始化第一类的子字典（300*7）
    temp=zeros(m1,size(cdat,2),k1);
    for i=1:k1
        dict=FDDL_INID(cdat,size(cdat,2),opts.wayInit);         % 100：每一类原子数
        temp(:,:,i)=dict;
    end
    Dict_ini      =    [Dict_ini temp];      %1：7列为第一类的初始化字典，8：14列为第二类的初始化子字典
    Dlabel_ini    =    [Dlabel_ini repmat(ci,[1 size(temp,2)])];
end


%    dict_path = ['E:\Matlab2016a\Aworkplace\2020classify_SR\Dictionary\D_ini' '.mat' ];
%    save(dict_path, 'Dict_ini');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 一类一类初始化编码系数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ini_par.tau         =     opts.lambda1;
ini_par.lambda      =     opts.lambda2;
ini_ipts.D          =     Dict_ini;            %初始化的大字典
 coef = zeros(size(Dict_ini,2),size(TrainDat,2),k1);
if size(Dict_ini,1)>size(Dict_ini,2)           %d = eigs(A,k) %返回k个最大特征值
      ini_par.c        =    1.05*eigs(Dict_ini(:,:,1)'*Dict_ini(:,:,1),1);
else
      ini_par.c        =    1.05*eigs(Dict_ini(:,:,1)*Dict_ini(:,:,1)',1);
end
for ci =  1:opts.nClass     %一类一类的初始化编码系数
    fprintf(['Initializing Coef:  Class ' num2str(ci) '\n']);
    ini_ipts.X      =    TrainDat(:,TrainLabel==ci,:);
    [ini_opts]      =    FDDL_INIC (ini_ipts,ini_par);
    coef(:,TrainLabel ==ci,:) =    ini_opts.A;   
    ert = ini_opts.ert;
end                                            
                                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Main loop of Dictionary Learning
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 Fish_par.dls        =     Dlabel_ini;   %字典标签
 Fish_ipts.D         =     Dict_ini;     %初始化字典
 Fish_ipts.trls      =     TrainLabel;    %样本类别标签
 Fish_par.tau        =     opts.lambda1;
 Fish_par.lambda2    =     opts.lambda2;
 Fish_nit            =     1;
 drls                =     Dlabel_ini;
 Fish_par.TrainDat   = TrainDat;
 CoefM=zeros(size(Fish_ipts.D,2),1,size(Fish_ipts.D,3));
 mu=0.0001;
 mu_1=0.5;
 % 张量逆，傅里叶变换求逆，在逆变换
 D=Fish_ipts.D;
 A=coef;
 
 V=zeros(size(A(TrainLabel ==ci,Fish_ipts.trls==ci,:)));
for ci=1:opts.nClass
    
    Zi=TrainDat(:,TrainLabel ==ci,:);
    Di=D(:,TrainLabel ==ci,:);
    Ai=A(TrainLabel ==ci,Fish_ipts.trls==ci,:);
    for iii=1:10
        DiT=Trans(Di);
        DDT=tprod(DiT,Di);         
%         D(:,TrainLabel ==ci,:)=1;
        E1=zeros(size(Di,2),size(Di,2),size(Di,3));
        E1(:,:,1)=eye(size(Di,2));
%         E1(:,:,2)=eye(size(Di,2));
%         E1(:,:,3)=eye(size(Di,2));
        B_1=inverse(DDT+mu_1.*E1);
        B_2=tprod(DiT,Zi)+mu_1.*Ai+V/2;
        B=tprod(B_1,B_2);
        Ai=soft(B-V/(2*mu_1),mu/(2*mu_1));
        AiT=Trans(Ai);
        iAAT=inverse(tprod(Ai,AiT));        
        Di=tprod(tprod(Zi,AiT),iAAT);

        V=V+0.1*mu_1.*(Ai-B);

    end
    Fish_ipts.D(:,TrainLabel ==ci,:)=Di;
    A(TrainLabel ==ci,Fish_ipts.trls==ci,:)=Ai;
end
    
    
    
    Dict = Fish_ipts.D; 
    Drls = drls;
%     dict_path = ['E:\Matlab2016a\Aworkplace\FDDL\Dictionary\bigD_' '.mat' ];
%     save(dict_path, 'Dict');
  
return;

%  -----------------------------------------------------------------------------------------------------------

for ci=1:opts.nClass
    
    Zi=TrainDat(:,TrainLabel ==ci,:);
    Di=D(:,TrainLabel ==ci,:);
    Ai=A(TrainLabel ==ci,Fish_ipts.trls==ci,:);
    for iii=1:10
        upadte_D();
        update_B1();
        upadte_Aii();        
        update_Lamda;
        update_B2();
        update_Aij();
        update_Lamda2;
        
    end
    Fish_ipts.D(:,TrainLabel ==ci,:)=Di;
    A(TrainLabel ==ci,Fish_ipts.trls==ci,:)=Ai;
end



