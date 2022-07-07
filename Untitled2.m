Xh = zeros(5^2,500,31);
Xl = zeros(4*5^2,500,31);

for j=1:31 % ¹âÆ×Ñ­»·
    XH=[];
    XL=[];
    for ii = 1:1
        temp = a(:,:,j);
        [H, L] = sample_patches(temp, 5, 500, 2);
        XH = [XH, H];
        XL = [XL, L];
    end
    Xh(:,:,j)=XH;
    Xl(:,:,j)=XL;
end


% [XX,label,cls_num] = patch_class(Xh, Xl,3,15);

