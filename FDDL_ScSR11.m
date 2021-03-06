function [hIm] = FDDL_ScSR11(lIm, up_scale, Dh, Dl, im_h,drls)
A = double(lIm);    
p = fspecial('gaussian', 5,0.8);                  
% p = fspecial('gaussian', 5,1);

p = p.^2;
p = p./sum(p(:)); 
im=convn(im_h, p, 'same');
[h,w,la] = size(im_h);
F=create_F();
mIm = imresize(A, up_scale, 'bicubic');
% t= hyperConvert3D((F*hyperConvert2D(im)), size(mIm,1), size(mIm,2));
[m1,n1,k1]=size(Dh);
[m2,n2,k2]=size(Dl);

[m3,n3,k3]=size(lIm);

norm_Dl = sqrt(sum(Dl.^2, 1)); 
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


s         = 2;  
b         = patch_size;
N         =  h-b+1;                        
M         =  w-b+1;
r         =  [1:s:N];           
r         =  [r r(end)+1:N];      
c         =  [1:s:M];
c         =  [c c(end)+1:M];
L         =  length(r)*length(c);   

X_m=fea_patch(b,lImfea_Z,r,c);
X=X_m;
X=normal(X);
X_m=is_nan(X_m);
load('dictionary_trained/center_patches_cluster_3.mat'); 
vec1  =vec;                             
set         =   1:size(X, 2);                              
L           =   size(set,2);                              
b2          =   size(X, 1);
X_t1=X(:,:,1);
for j = 1 : L        
%     jj = set(j);
%     vv= vec(:,1);    
%     yy= X(1,set(j));   
    dis   =   (vec(:, 1) -  X_t1(1, set(j))).^2;  
    for i = 2 : b2
        dis  =  dis + (vec(:, i)-X_t1(i, set(j))).^2;   
    end
    [~,ind]      =   min( dis );  
    cls_idx( set(j) )   =   ind;  
end
cls_idx = cls_idx';
[s_idx, seg]    =   Proc_cls_idx( cls_idx ); 

Y = zeros( b*b, L , k1);
Fish_ipts.Dl        = Dl;      
Fish_par.dls        = drls;  
Fish_par.tau        = 0.005; 
Fish_par.cls_num    = length(unique(drls));   
Fish_par.nIter      = 100;
% Fish_par.c = 1.05*eigs(Dl(:,:,1)*Dl(:,:,1)',1);      
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
        D_n = Dl(:,drls~=ci,:); 
        Ain = A(drls~=ci,idx ,:);
        V=zeros(size(A(drls==ci,idx ,:)));
        VV1=zeros(size(A(drls~=ci,idx ,:)));
        Aii=A(drls==ci,idx ,:);
        Zi=X(:,idx,:);
        for iii=1:2
            
            Zi=is_nan(Zi);
    %         Di=D(:,TrainLabel ==ci,:);
            DnT = Trans(D_n);
            B1_1 = tprod(DnT, D_n);
            EEE = zeros(size(B1_1));
            EEE(:,:,1)=eye(size(B1_1,1));
            B1_2 = tprod(DnT,Zi-tprod(Di, Aii))*mu + mu_1*(Ain+mu*VV1/(2*mu_1));
            B1 = tprod(inverse(2*B1_1+mu_1*EEE),B1_2);
            Ain = soft(B1-mu*VV1/(2*1), 0.05);
            VV1=VV1+mu_1*(Ain-B1);
            DiT=Trans(Di);
            DDT=tprod(DiT,Di);
            
    %         D(:,TrainLabel ==ci,:)=1;
            E1=zeros(size(Di,2),size(Di,2),size(Di,3));
            E1(:,:,1)=eye(size(Di,2));
            B_1=inverse(DDT+mu_1.*E1);
            B_2=tprod(DiT,Zi)+mu_1.*Aii+V/2;
            B=tprod(B_1,B_2);
            Aii=soft(B-V/(2*mu_1),mu);
            Aii=is_nan(Aii);
            AiT=Trans(Aii);            
            iAAT=inverse(tprod(Aii,AiT));
            Di=tprod(tprod(Zi,AiT),iAAT);
            V=V+mu_1.*(Aii-B);
        end
        A(drls==ci,idx,:)=Aii;
        A(drls~=ci,idx,:)=Ain;
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

im_out   =  zeros(h,w,spec);           
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
him=im1;

% hIm=Log_prox_tnn( him, eta/2/mu,up_scale ,t);
hIm=Log_prox_tnn( him, eta/2/mu);


