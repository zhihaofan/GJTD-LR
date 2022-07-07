function [hIm]=fidelity(im_l,im_h,maxIter,upscale)
% k1=size(im_l,3);
% [row_l, col_l, spec_l] = size(im_l);
% [row_h, col_h, spec_h] = size(im_h);

p = fspecial('gaussian', 5, 1.56);
p = p.^2;
p = p./sum(p(:)); 

% p_3d=repmat(p,[1,1,k1]);


im_l = double(im_l);
im_h = double(im_h);

for ii = 1:maxIter,
    im_l_s = imresize(im_h,1/upscale,'bicubic');
    im_diff = im_l - im_l_s;
    
    im_diff = imresize(im_diff, upscale, 'bicubic');
    im_h = im_h + convn(im_diff, p, 'same');
end
hIm=im_h;

