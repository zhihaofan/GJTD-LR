function [XH, XL] = t_patch_pruning(Xh, Xl, threshold)
[a1,a2,a3]=size(Xh);
[b1,b2,b3]=size(Xl);

pvars = var(Xh(:,:,1), 0, 1);

idx = pvars > threshold;
I=find(idx(:)~=0);
num=length(I) ;
XH=zeros(a1,num,a3);
XL=zeros(b1,num,b3);
for i=1:a3
    tem=Xh(:,:,i);
    XH(:,:,i)=tem(:,idx);
    tem=Xl(:,:,i);
    XL(:,:,i)=tem(:,idx);
end