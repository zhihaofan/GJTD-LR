function [Dict,Drls,coef] = FDDL(TrainDat,TrainLabel,opts)
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
Dict_ini  =  [];  
Dlabel_ini = []; 

for ci = 1:opts.nClass
    cdat          =    TrainDat(:,TrainLabel==ci,:);   
%     dict          =    FDDL_INID(cdat,size(cdat,2),opts.wayInit);  
    temp=zeros(m1,size(cdat,2),k1);
    for i=1:k1
        dict=FDDL_INID(cdat,size(cdat,2),opts.wayInit);         
        temp(:,:,i)=dict;
    end
    Dict_ini      =    [Dict_ini temp];     
    Dlabel_ini    =    [Dlabel_ini repmat(ci,[1 size(temp,2)])];
end


%    dict_path = ['E:\Matlab2016a\Aworkplace\2020classify_SR\Dictionary\D_ini' '.mat' ];
%    save(dict_path, 'Dict_ini');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ini_par.tau         =     opts.lambda1;
ini_par.lambda      =     opts.lambda2;
ini_ipts.D          =     Dict_ini;            
 coef = zeros(size(Dict_ini,2),size(TrainDat,2),k1);
if size(Dict_ini,1)>size(Dict_ini,2)          
      ini_par.c        =    1.05*eigs(Dict_ini(:,:,1)'*Dict_ini(:,:,1),1);
else
      ini_par.c        =    1.05*eigs(Dict_ini(:,:,1)*Dict_ini(:,:,1)',1);
end
for ci =  1:opts.nClass    
    fprintf(['Initializing Coef:  Class ' num2str(ci) '\n']);
    ini_ipts.X      =    TrainDat(:,TrainLabel==ci,:);
    [ini_opts]      =    FDDL_INIC (ini_ipts,ini_par);
    coef(:,TrainLabel ==ci,:) =    ini_opts.A;   
    ert = ini_opts.ert;
end                                            
                                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Main loop of Dictionary Learning
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 Fish_par.dls        =     Dlabel_ini;  
 Fish_ipts.D         =     Dict_ini;   
 Fish_ipts.trls      =     TrainLabel;   
 Fish_par.tau        =     opts.lambda1;
 Fish_par.lambda2    =     opts.lambda2;
 Fish_nit            =     1;
 drls                =     Dlabel_ini;
 Fish_par.TrainDat   = TrainDat;
 CoefM=zeros(size(Fish_ipts.D,2),1,size(Fish_ipts.D,3));
 mu=0.0001;
 mu_1=0.5;
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



