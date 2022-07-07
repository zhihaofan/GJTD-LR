%clustering 最终返回聚类的数目（cls_num）每个块块所属类别（类别序号cls_idx),类中心（vec)
function   [cls_idx,vec,cls_num]  =  Clustering(Y, cls_num, itn)
Y         =   Y';                %  一行为一块。本来是25*1000，转置为1000*25，一行即为一块
[L b2]    =   size(Y);           % 总行数，块的数量
P         =   randperm(L);       %将图像块序号随机打乱 即：产生1~图像块数目间的随机顺序的序号
P2        =   P(1:cls_num);      %随机选取3个图像的序号
vec       =   Y(P2(1:end), :);   %10*100 一行为一个块；选取10个图像块序号对应的图像块，作为分别作为10个类的中心
m_num     =   2;               %控制每个类中的块数目大于
%根据类中心将其他图像块分类
for i = 1 : itn     %k-means 迭代15次
    cnt       =  zeros(1, cls_num);    
    v_dis    =   zeros(L, cls_num);
    for  k = 1 : cls_num     %k = 1 所有图像块到类中心1的距离; 
        v_dis(:, k) = (Y(:,1) - vec(k,1)).^2;
        for c = 2:b2
            v_dis(:,k) =  v_dis(:,k) + (Y(:,c) - vec(k,c)).^2; 
        end              %块数*10（类别数）                      %第一行：第一块到所有类中心的距离
    end                                                          %第二行：第二块到所有类类中心的距离
    
    [val,cls_idx]     =   min(v_dis, [], 2);    %cls_idx 图像块距离哪个类别距离最小，就保存该类别的序号。求（A,[],2）所组成数组的最小值, [t,C]=min(A,[],2)该函数返回最小值到t中,返回最小值的序号到C中
   
    [s_idx, seg]   =  Proc_cls_idx( cls_idx );  %s_idx类别序号（64个）排序后对应的图像块序号(74047*1)
    for  k  =  1 : length(seg)-1                %seg每各类最后一个图像块的为第多少个图像块，找到类最后一块的序号，以便于分类
        idx    =   s_idx(seg(k)+1:seg(k+1));    %idx:k=1时，存放属于第一类的图像块共1271个（1~1271元素序号对应图像块的序号）
        cls    =   cls_idx(idx(1));             %类标签
        vec(cls,:)    =   mean(Y(idx, :));      %对每个类别求均值，获得新的类中心（像素值？）第k 类的新的类中心
        cnt(cls)      =   length(idx);          %每一个类别的含有的图像块数目
    end
    
                         %一直循环1：11次，对这些图像块进行聚类（改中心，算距离，聚类）
    if (i==itn-2)        %i= 13时    
        [val ind]  =  min( cnt );  %找出这64个类中，包含块数最少的块 
%         while (val<m_num) && (cls_num>=40)  %while 循环的作用，如果分类结果中，类中块的数目<100那么就将这样的类别省去
          while (val<m_num)
            vec(ind, :)    =  [];
            cls_num       =  cls_num - 1;   
            cnt(ind)      =  [];
            [val  ind]    =  min(cnt);
           end        
    end
end
end