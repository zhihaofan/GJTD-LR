function [Dict,Drls,coef] = new_FDDL_SR(TrainDat,TrainLabel,opts)
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
% for ci =  1:opts.nClass    
%     fprintf(['Initializing Coef:  Class ' num2str(ci) '\n']);
%     ini_ipts.X      =    TrainDat(:,TrainLabel==ci,:);
%     [ini_opts]      =    FDDL_INIC (ini_ipts,ini_par);
%     coef(:,TrainLabel ==ci,:) =    ini_opts.A;   
%     ert = ini_opts.ert;
% end
load('./init.mat');
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
 
par_set.mu_1 = 0.0001;
par_set.mu_2 = 0.5;
size0 = size(A(TrainLabel == ci,Fish_ipts.trls == ci, :));
Init_par.Lamada1 = zeros(opts.nClass, size0(1), size0(2), size0(3));
size1 = size(A(TrainLabel ~= ci,Fish_ipts.trls == ci, :));
Init_par.Lamada2 = zeros(opts.nClass, size1(1), size1(2), size1(3));
iii = 0;
while iii < 30
    for ci=1:opts.nClass

        Zi=TrainDat(:,TrainLabel ==ci,:);
        Di=Fish_ipts.D(:,TrainLabel ==ci,:);
        Ai=A(Fish_par.dls ==ci,Fish_ipts.trls==ci,:);
        D = update_D(par_set, D, Zi, A, ci, TrainLabel, Fish_ipts);
        EI = zeros(size(Di, 2), size(Di, 2), size(Di, 3));
        EI(:, :, 1) = eye(size(Di, 2));
        Aija = A(TrainLabel ~= ci,Fish_ipts.trls == ci,:);
        Dja = D(:,TrainLabel ~=ci,:);
        H = Zi - tprod(Dja, Aija);
        Aii = A(TrainLabel ==ci,Fish_ipts.trls==ci,:);
        B1 = update_B1(par_set, D, EI, H, Zi, Aii, Init_par.Lamada1, ci, TrainLabel);
        Aii = update_Aii(par_set, B1, Init_par.Lamada1, ci);
        A(TrainLabel ==ci,Fish_ipts.trls==ci,:) = Aii;
        Init_par.Lamada1(ci, :, :, :) = Init_par.Lamada1(ci, :, :, :) + reshape(par_set.mu_1 * (Aii - B1), 1, size0(1), size0(2), size0(3));
        EI1 = zeros(size(Dja, 2), size(Dja, 2), size(Dja, 3));
        EI1(:, :, 1) = eye(size(Dja, 2));
        H1 = Zi - tprod(Dja,A(TrainLabel ~=ci,Fish_ipts.trls==ci,:));
        B2 = update_B2(par_set, Dja, Aija, EI1, H1, Init_par.Lamada2, ci);
        Aija = update_Aij(par_set, B2, Init_par.Lamada2, ci);
        A(TrainLabel ~=ci,Fish_ipts.trls==ci,:) = Aija;
        Init_par.Lamada2(ci, :, :, :) = Init_par.Lamada2(ci, :, :, :) + reshape(par_set.mu_1 * (Aija - B2), 1, size(Aija, 1), size(Aija, 2), size(Aija, 3));
%         cnt = 0;
%         for ij = 1:opts.nClass
%             if ij == ci
%                 continue
%             end
%             cnt = cnt +1;
%             Aij = A(TrainLabel ==ij,Fish_ipts.trls==ij,:);
%             Lamada_t = reshape(Init_par.Lamada2(ci, cnt, :, :, :), size(Init_par.Lamada1));
%             Lamada_t = Lamada_t + par_set.mu_1 * (Aij - B2);
%             Init_par.Lamada2(ci, cnt, :, :, :) = reshape(Lamada_t, [1, 1, size(Init_par.Lamada1)]);
%         end

        Fish_ipts.D=D;
%         A(TrainLabel ==ci,Fish_ipts.trls==ci,:)=Ai;
    end
    iii = iii + 1;
end
    
    
Dict = Fish_ipts.D; 
Drls = drls;
%     dict_path = ['E:\Matlab2016a\Aworkplace\FDDL\Dictionary\bigD_' '.mat' ];
%     save(dict_path, 'Dict');
  
return;

%  -----------------------------------------------------------------------------------------------------------




