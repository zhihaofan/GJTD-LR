function  [s_idx, seg]    =  Proc_cls_idx( cls_idx )    %类的序号  属于的类别 cls_idx 每个块所属于的类别
[idx  s_idx]    =  sort(cls_idx);          %idx 存放类序号排序结果第一类，第二类...，
                                           %s_idx 排序结果索引块的位置，属于第一类得块，属于第二类的块，
idx2   =  idx(1:end-1) - idx(2:end);   %错位相减
seq    =  find(idx2);    %63*1 错位相减后有263个类序号差非零，返回对应块的索引位置

seg    =  [0; seq; length(cls_idx)];   