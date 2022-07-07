%--------------------------------------------------------------------------
function [pos_arr, wei_arr]  =  Block_matching(im, par)
S         =  par.s1;   %规定搜索范围的大小
f         =  par.win;
f2        =  f^2;
nv        =  par.nblk;
s         =  par.step;     %为保证6*6的图像块重叠4个像素
hp        =  max(12*par.nSig, par.hp);  %hp

N         =  size(im,1)-f+1;
M         =  size(im,2)-f+1;
r         =  [1:s:N];
r         =  [r r(end)+1:N];
c         =  [1:s:M];
c         =  [c c(end)+1:M];
L         =  N*M; %共64009块
X         =  zeros(f*f, L);

%获取6*6叠5个像的图像块结果：获得64009=253*253个6*6的图像块存放于X当中
k    =  0;
for i  = 1:f
    for j  = 1:f
        k    =  k+1;
        blk  =  im(i:end-f+i,j:end-f+j);
        X(k,:) =  blk(:)';
    end
end


I     =   (1:L);
I     =   reshape(I, N, M);
N1    =   length(r);
M1    =   length(c);   %6*6的块叠4个像素分块，获得的列数127
pos_arr   =  zeros(nv, N1*M1 );  %在6*6叠4分块基础上，用于存放12个最相似的图像块，
wei_arr   =  zeros(nv, N1*M1 );  %每个相似块对应的权重
X         =  X';                 %64009*36 一行为一个块，分块方式：6*6叠5个像素
%
for  i  =  1 : N1
    for  j  =  1 : M1
        row     =   r(i);      %r中存放横着（行上）取1，3，5，7，，，127图像像素的序号，取奇数序号的像素
        col     =   c(j);      %竖着（列上）
        off     =  (col-1)*N + row;  %在6*6叠4分割的基础上，off:图像像素的序号，确定取哪些像素：1，507，、、
        off1    =  (j-1)*N1 + i;                           %off1:分割图像块的序号，最大为16129（127*127）块按6*6叠4 分，共得到16129个块，
        rmin    =   max( row-S, 1 );    %确定搜索区域的图像像素序号
        rmax    =   min( row+S, N );
        cmin    =   max( col-S, 1 );
        cmax    =   min( col+S, M );
        idx     =   I(rmin:rmax, cmin:cmax);
        idx     =   idx(:);     %确定的搜索范围是像素的序号，对应找每个像素对应块，B
        B       =   X(idx, :);     %取搜索框对应像素序列号的对应块的对应像素值；64009*36 一行为一个块，分块方式：6*6叠5个像素   
        v       =   X(off, :);     %由图像像素的序号，找到对应的6*6的图像块
        dis     =   (B(:,1) - v(1)).^2;  %计算像素序号对应图像块与搜索范围内块的距离
        for k = 2:f2
            dis   =  dis + (B(:,k) - v(k)).^2;
        end
        dis   =  dis./f2;
        [val,ind]   =  sort(dis);    %距离的值，ind: 排序后，距离对应的图像块（1~676 的）排序序号   
        dis(ind(1))  =  dis(ind(2));  %第一个距离设置为除它自己本身外，距离它最近的图像块的距离     
        wei         =  exp( -dis(ind(1:nv))./hp );  %权重矩阵12*1
        wei         =  wei./(sum(wei)+eps);         %权重归一化操作
        indc        =  idx( ind(1:nv) );            %取出距离最近的12块图像块序号
        pos_arr(:,off1)  =  indc;               %第一次循环时；pos_arr的第一列存放最相似图像块的序号，第二次循环时，第128列存放前十二个相似块 序号
        wei_arr(:,off1)  =  wei;               %第一次循环时；wei_arr的第一列存放最相似图像块的权重
    end   
end
pos_arr  =  pos_arr';
wei_arr  =  wei_arr';
