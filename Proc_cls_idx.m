function  [s_idx, seg]    =  Proc_cls_idx( cls_idx )    %������  ���ڵ���� cls_idx ÿ���������ڵ����
[idx  s_idx]    =  sort(cls_idx);          %idx ����������������һ�࣬�ڶ���...��
                                           %s_idx �������������λ�ã����ڵ�һ��ÿ飬���ڵڶ���Ŀ飬
idx2   =  idx(1:end-1) - idx(2:end);   %��λ���
seq    =  find(idx2);    %63*1 ��λ�������263������Ų���㣬���ض�Ӧ�������λ��

seg    =  [0; seq; length(cls_idx)];   