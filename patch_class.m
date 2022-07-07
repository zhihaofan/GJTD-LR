function [XX,label,cls_num] = patch_class(Xh, Xl,class_num,k_iter)
addpath(genpath('RegularizedSC'));
hDim = size(Xh, 1);
lDim = size(Xl, 1);
spec = size(Xh,3);
for i=1:spec
    Xh(:,:,i) = Xh(:,:,i)./repmat(sqrt(sum(Xh(:,:,i).^2)), size(Xh(:,:,i), 1), 1);     %���зֱ�Ըߵͷֱ���������淶������һ��
    Xl(:,:,i) = Xl(:,:,i)./repmat(sqrt(sum(Xl(:,:,i).^2)), size(Xl(:,:,i), 1), 1);
end
% joint learning of the dictionary
xh=sqrt(hDim)*Xh;                     %���淶���������е���������5��
xl=sqrt(lDim)*Xl;

                                                              
%% ��������ಢȡ1000��
 [cls_idx,vec,cls_num]  =  Clustering(xl(:,:,1), class_num, k_iter);   %cls_idx �����ÿ�����������ڵ����k-means���࣬�����࣬����15��
 %��������vec 25*10/9/8�������� ���������������ؽ�ʱ��Ƚ�������������һ��
% vec_path = ['dictionary_trained\center_patches' num2str(cls_num) '.mat' ];
save('./vec.mat ', 'vec');  
[s_idx, seg]   =  Proc_cls_idx( cls_idx );  %s_idx�����ţ�����������Ӧ��ͼ�����ţ���һ�࣬�ڶ��࣬����ͼ������      seg��ǰ�����һ��������ֵ
Xh_new = []; 
Xl_new = []; 
label = [];
 
 % �Էֿ��X���з��࣬���������
for  i  =  1 : length(seg)-1                %ȡ��ÿ���Ӧ�Ŀ�
    idx    =   s_idx(seg(i)+1:seg(i+1));    %i=1ʱ��ȡ����һ���ͼ����Ӧ��ͼ������
    idx1   = idx(1:200,:);                  %ÿ������ȡͼ�����Ŀ 200����
    xl_l    = xl(:, idx1,:);                   %i=1 ʱ��ȡ����Ӧ�ĵͣ��飩��Ϊ��1��������
    xh_h    = xh(:,idx1,:);                    %��Ӧ��ȡ���߿�
                                      
    Xh_new = [Xh_new,xh_h];
    Xl_new = [Xl_new,xl_l];
    clss = cls_idx(idx1);                   %ȡ����Ӧ������ǩi
    clss = clss';
    label = [label clss];      %����ͼ�������ǩ 
    
end  

XX = [Xh_new; Xl_new];                  %�������е�����ѵ������

% for i=1:spec
%     XX(:,:,i) = XX(:,:,i)./repmat(sqrt(sum(XX(:,:,i).^2, 1)), hDim+lDim, 1);    % ��һ������ֵ�ѵ������125*1000��1000��������
% end
