function [hIm] = new_FDDL_ScSR_par(lIm, up_scale, Dh, Dl, im_h,drls, par_1)
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
% load('dictionary_trained/center_patches_cluster_3.mat'); 
load('vec_3.mat'); 
vec1  =vec;                           
set         =   1:size(X, 2);                              
L           =   size(set,2);                               
b2          =   size(X, 1);
X_t1=X(:,:,1);
for j = 1 : L       
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
%         Fish_par.tau        = 0.0015; 
Fish_par.cls_num    = length(unique(drls));  
Fish_par.nIter      = 100;
% Fish_par.c = 1.05*eigs(Dl(:,:,1)*Dl(:,:,1)',1);    
Fish_nit=1;
eta=0.0025;
mu=0.0025;
mu_1=0.05;

% Initialize
Init_par.Lamada1 = cell(length(seg)-1, 1);
for ci=1:length(seg)-1
    idx = s_idx(seg(ci)+1:seg(ci+1));
    Init_par.Lamada1{ci} = zeros(length(Fish_par.dls == ci), length(idx), k3);
end


A_l = preprocess_A(X, Dl, seg, Fish_par, s_idx);
while Fish_nit <= 5
    if Fish_nit == 1
        A_iter = A_l;
    end
    [A_iter, Init_par] = update_SR_A_par(Init_par, Y, Dh, A_iter, Fish_par, seg, s_idx, A_l, par_1);
    for ci=1:length(seg)-1  % class
        
        idx    =   s_idx(seg(ci)+1:seg(ci+1));
        
        Mi = merge_m(Dh, A_iter, ci, idx, Fish_par);

        Xi = update_SR_Xi(Mi);
        Y(:,idx,:)=Xi;
    end
    Fish_nit = Fish_nit + 1;
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
hIm=im1;

% hIm=Log_prox_tnn( him, eta/2/mu,up_scale ,t);
% hIm=Log_prox_tnn( him, eta/2/mu);




return;
