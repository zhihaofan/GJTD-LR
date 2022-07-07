function [hIm] = FDDL_ScSR(lIm, up_scale, Dh, Dl, im_h,drls)
%% 迭代次数控制：Fish_par.nIter
%% 进度：109行
A = double(lIm);     %待重建图像
p = fspecial('gaussian', 5,0.8);                  
% p = fspecial('gaussian', 5,1);

p = p.^2;
p = p./sum(p(:)); 
im=convn(im_h, p, 'same');
[h,w,la] = size(im_h);
F=create_F();
%将图像分块并提取图像块特征
mIm = imresize(A, up_scale, 'bicubic');
t= hyperConvert3D((F*hyperConvert2D(im)), size(mIm,1), size(mIm,2));
[m1,n1,k1]=size(Dh);
[m2,n2,k2]=size(Dl);

[m3,n3,k3]=size(lIm);

norm_Dl = sqrt(sum(Dl.^2, 1)) + 1e-15;
Dl = Dl./repmat(norm_Dl, size(Dl, 1), 1);

patch_size = sqrt(size(Dh, 1));

% bicubic interpolation of the low-resolution image
% mIm = double(imresize(lIm, up_scale, 'bicubic')); 

hIm = zeros(size(mIm));
cntMat = zeros(size(mIm));

[h, w, spec] = size(mIm);

% extract low-resolution image features\
lImfea_Z=zeros(up_scale*m3,up_scale*n3,4,k1);
for i=1:k1
    lImfea_Z(:,:,:,i) = extr_lIm_fea(mIm(:,:,i));
end


s         = 2;     %5*5 叠3
b         = patch_size;
N         =  h-b+1;                          %取得的最后一个像素为253，确保最后一个块为6*6的大小
M         =  w-b+1;
r         =  [1:s:N];              %要取的，图像行像素的序号
r         =  [r r(end)+1:N];       % 1~253  把循环生成的数值放到数组r 中，为后面循环获取图像块做准备。防止取块时，少一部分,补【r(end)+1:N】
c         =  [1:s:M];
c         =  [c c(end)+1:M];
L         =  length(r)*length(c);   %253*253=64009 图像像素序号的总的长度  当s=2时，L=127*127=16129块间重叠2

X_m=fea_patch(b,lImfea_Z,r,c);
X=X_m;
X=normal(X);
X_m=is_nan(X_m);
%% 聚类标签
load('dictionary_trained/center_patches_cluster_3.mat'); 
% load('./vec.mat');
vec1  =vec;                             %一行是一个类别中心
set         =   1:size(X, 2);                                %重建图像块个数
L           =   size(set,2);                                 %重建图像块个数    
b2          =   size(X, 1);
X_t1=X(:,:,1);
for j = 1 : L         %找到L个块分别属于哪一类
%     jj = set(j);
%     vv= vec(:,1);      %j = 1,得出第一个块到10个类中心的距离
%     yy= X(1,set(j));   %J= 2，得出第二个块到10个类中心的距离
    dis   =   (vec(:, 1) -  X_t1(1, set(j))).^2;    %set（1*10124）中存放1~10124个图像块的序号
    for i = 2 : b2
        dis  =  dis + (vec(:, i)-X_t1(i, set(j))).^2;   %i=1时，61*1
    end
    [~,ind]      =   min( dis );  %%ind 类别序号  set(j)：图像块序号
    cls_idx( set(j) )   =   ind;    %cls_idx中存放1~16129个块所属于的类别
end
cls_idx = cls_idx';
[s_idx, seg]    =   Proc_cls_idx( cls_idx );  %1~16129（127*127）个块所属类别（类别序号cls_idx)

Y = zeros( b*b, L , k1);
Fish_ipts.Dl        = Dl;      %100*1000
Fish_par.dls        = drls;   %字典类标签
Fish_par.tau        = 0.005;  %稀疏正则化参数
%         Fish_par.tau        = 0.0015;  %稀疏正则化参数
Fish_par.cls_num    = length(unique(drls));    %子字典类别数
Fish_par.nIter      = 100;
% Fish_par.c = 1.05*eigs(Dl(:,:,1)*Dl(:,:,1)',1);      % 参数sigma
Fish_nit=1;
eta=0.0025;
mu=0.0025;
mu_1=0.05;


A=zeros(size(Dl,2),size(X,2),k1);
while Fish_nit <= 1
    for ci=1:length(seg)-1  % class
        
        idx    =   s_idx(seg(ci)+1:seg(ci+1));
        Di=Dl(:,drls==ci,:); 
        Di=is_nan(Di);
        V=zeros(size(A(drls==ci,idx ,:)));
        Zi=X(:,idx,:);
        Aii=A(drls==ci,idx ,:);
        for iii=1:2
            
            Zi=is_nan(Zi);
    %         Di=D(:,TrainLabel ==ci,:);
            
            DiT=Trans(Di);
            DDT=tprod(DiT,Di);
            
    %         D(:,TrainLabel ==ci,:)=1;
            E1=zeros(size(Di,2),size(Di,2),size(Di,3));
            E1(:,:,1)=eye(size(Di,2));
            B_1=inverse(DDT+mu_1.*E1);
            B_2=tprod(DiT,Zi)+mu_1.*Aii+V/2;
            B=tprod(B_1,B_2);
%             Aii=soft(B-V/(2*mu_1),mu^2/(2*mu_1));
            Aii=soft(B-V/(2*mu_1),mu);
            Aii=is_nan(Aii);
            AiT=Trans(Aii);            
            iAAT=inverse(tprod(Aii,AiT));
            Di=tprod(tprod(Zi,AiT),iAAT);
            V=V+mu_1.*(Aii-B);
        end
        A(drls==ci,idx,:)=Aii;
        Ai=A(:,idx,:);
        Temp=tprod(Dh,Ai);
%         [a1,a2,a3]=size(Temp);
%         Temp_s=reshape(Temp,[a1*a2 a3]);
%         V22=WNNM( Temp_s, D2/2/mu, nsig2 );
%         V22_s=reshape(V22,a1,a2,a3);
%         V22_s=Log_prox_tnn( Temp, eta/2/mu );
%         Y(:,idx,:)=V22_s;
        Y(:,idx,:)=Temp;
    end
    Fish_nit=Fish_nit+1;
end

im_out   =  zeros(h,w,spec);            %存放重建图像 256*256
im_wei   =  zeros(h,w,spec);
k        =  0;
N        = length(r);
M        = length(c);
num=0;
for ii=1:N
    for jj=1:M
        num=num+1;
        kk= reshape(Y(:,num,:),[patch_size patch_size k1]);
        xx=r(ii);
        yy=c(jj);
        temp=mIm(yy:yy+patch_size-1,xx:xx+patch_size-1,:);
        mMean=mean(temp(:));
        im_out(yy:yy+patch_size-1,xx:xx+patch_size-1,:)=im_out(yy:yy+patch_size-1,xx:xx+patch_size-1,:)+kk+mMean;
        im_wei(yy:yy+patch_size-1,xx:xx+patch_size-1,:)=im_wei(yy:yy+patch_size-1,xx:xx+patch_size-1,:)+1;
    end
end


im1    =  im_out./(im_wei+eps);
% him=(2*im1 + 1*im_h)/3;
him=im1;

% hIm=Log_prox_tnn( him, eta/2/mu,up_scale ,t);
hIm=Log_prox_tnn( him, eta/2/mu);


