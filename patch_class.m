function [XX,label,cls_num] = patch_class(Xh, Xl,class_num,k_iter)
addpath(genpath('RegularizedSC'));
hDim = size(Xh, 1);
lDim = size(Xl, 1);
spec = size(Xh,3);
for i=1:spec
    Xh(:,:,i) = Xh(:,:,i)./repmat(sqrt(sum(Xh(:,:,i).^2)), size(Xh(:,:,i), 1), 1);     %按列分别对高低分辨率样本块规范化，归一化
    Xl(:,:,i) = Xl(:,:,i)./repmat(sqrt(sum(Xl(:,:,i).^2)), size(Xl(:,:,i), 1), 1);
end
% joint learning of the dictionary
xh=sqrt(hDim)*Xh;                     %将规范化的样本中的数扩大了5倍
xl=sqrt(lDim)*Xl;

                                                              
%% 样本块分类并取1000块
 [cls_idx,vec,cls_num]  =  Clustering(xl(:,:,1), class_num, k_iter);   %cls_idx 聚类后每个块所以属于的类别、k-means聚类，分三类，迭代15次
 %将类中心vec 25*10/9/8个类中心 保存起来，方便重建时候比较样本块属于那一类
% vec_path = ['dictionary_trained\center_patches' num2str(cls_num) '.mat' ];
save('./vec.mat ', 'vec');  
[s_idx, seg]   =  Proc_cls_idx( cls_idx );  %s_idx类别序号（个）排序后对应的图像块序号；第一类，第二类，。。图像块序号      seg：前后类别不一样的索引值
Xh_new = []; 
Xl_new = []; 
label = [];
 
 % 对分块的X进行分类，按类别排序
for  i  =  1 : length(seg)-1                %取出每类对应的块
    idx    =   s_idx(seg(i)+1:seg(i+1));    %i=1时，取出第一类的图像块对应的图像块序号
    idx1   = idx(1:200,:);                  %每个类所取图像块数目 200个块
    xl_l    = xl(:, idx1,:);                   %i=1 时，取出对应的低（块）作为第1类样本块
    xh_h    = xh(:,idx1,:);                    %对应低取出高块
                                      
    Xh_new = [Xh_new,xh_h];
    Xl_new = [Xl_new,xl_l];
    clss = cls_idx(idx1);                   %取出对应的类别标签i
    clss = clss';
    label = [label clss];      %所有图像块的类标签 
    
end  

XX = [Xh_new; Xl_new];                  %按类排列的联合训练样本

% for i=1:spec
%     XX(:,:,i) = XX(:,:,i)./repmat(sqrt(sum(XX(:,:,i).^2, 1)), hDim+lDim, 1);    % 归一化后的字典训练样本125*1000（1000个样本）
% end
