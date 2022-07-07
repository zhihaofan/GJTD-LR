function [XX,label,cls_num] = patch_class(Xh, Xl,class_num,k_iter)
addpath(genpath('RegularizedSC'));
hDim = size(Xh, 1);
lDim = size(Xl, 1);
spec = size(Xh,3);
for i=1:spec
    Xh(:,:,i) = Xh(:,:,i)./repmat(sqrt(sum(Xh(:,:,i).^2)), size(Xh(:,:,i), 1), 1);    
    Xl(:,:,i) = Xl(:,:,i)./repmat(sqrt(sum(Xl(:,:,i).^2)), size(Xl(:,:,i), 1), 1);
end
% joint learning of the dictionary
xh=sqrt(hDim)*Xh;                   
xl=sqrt(lDim)*Xl;

                                                              
 [cls_idx,vec,cls_num]  =  Clustering(xl(:,:,1), class_num, k_iter);   
% vec_path = ['dictionary_trained\center_patches' num2str(cls_num) '.mat' ];
save('./vec.mat ', 'vec');  
[s_idx, seg]   =  Proc_cls_idx( cls_idx ); 
Xh_new = []; 
Xl_new = []; 
label = [];
 
for  i  =  1 : length(seg)-1             
    idx    =   s_idx(seg(i)+1:seg(i+1)); 
    idx1   = idx(1:200,:);               
    xl_l    = xl(:, idx1,:);                 
    xh_h    = xh(:,idx1,:);               
                                      
    Xh_new = [Xh_new,xh_h];
    Xl_new = [Xl_new,xl_l];
    clss = cls_idx(idx1);                  
    clss = clss';
    label = [label clss];      
    
end  

XX = [Xh_new; Xl_new];                 

% for i=1:spec
%     XX(:,:,i) = XX(:,:,i)./repmat(sqrt(sum(XX(:,:,i).^2, 1)), hDim+lDim, 1);    
% end
