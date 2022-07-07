% =========================================================================
clear all; clc;     
tic;
addpath(genpath('quality')); 
addpath(genpath('LT'));

% load('Data/data.mat');
load('Data/face_ms.mat');

a=imresize(a,0.5,'bicubic');
% t=imresize(a,4,'bicubic');
% save('Data/aa','a')
% aa=aa.aa;

a=a(8:247,8:247,:);

% a=guiyihua(a,255,0);
% im_h=aa;
im_h=a;
% im_h=imresize(aa,0.5,'bicubic');
% set parameters
lambda = 0.1;                   
overlap = 3;                   
up_scale = 2;                   
maxIter = 50;                   
patch_size = 5;
% load dictionary
%  load('inidictionary&center&patches/s1_low1.2_D_last1000_0.15_5.mat');

tic;
tem=imresize(im_h,1/up_scale,'bicubic');
% aaaa = imresize(tem,up_scale,'bicubic');



p = fspecial('gaussian', 7,2);
p = p.^2;
p = p./sum(p(:)); 
im_l=convn(tem, p, 'same');
im_bic=imresize(im_l,up_scale,'bicubic');

% im_l=imresize(tem,2,'bicubic');
% im_l =imresize(im_h,1/up_scale,'bicubic');
toc;


% % load('Dictionary/3D_600_dim_31_class3.mat');
load('./new_Dict_10iter.mat');

D = Dh;
drls = zeros(1,size(D,2));
    s = 0;
for i =1:200:size(Dl,2);    
    s = s+1;
    drls(:,i:i+199) = s;
end

im = zeros(size(im_h));
[h, w, c] = size(tem);
im_patch = 128;
step = im_patch / up_scale;
stride = 2;
for x=0:(step - stride):h-1
    if x + step > h
        x_begin = h - step;
    else
        x_begin = x;
    end
    for y=0:(step - stride):w-1
        if y + step > w
            y_begin = h - step;
        else
            y_begin = y;
        end
        tem_patch = tem(x_begin + 1:x_begin+step, y_begin + 1:y_begin+step,:);
        im_h_patch = im(x_begin*up_scale + 1:x_begin*up_scale + im_patch, y_begin*up_scale+1:y_begin*up_scale + im_patch, :);
        im(x_begin*up_scale + 1:x_begin*up_scale + im_patch, y_begin*up_scale+1:y_begin*up_scale + im_patch, :) = new_FDDL_ScSR(tem_patch, up_scale, Dh, Dl, im_h_patch, drls);
    end
end
% im = new_FDDL_ScSR(tem, up_scale, Dh, Dl, im_h, drls);



% im=hIm;
toc;
im=is_nan(im);


% [image]=fidelity(im_l_ycbcr,im_h_ycbcr,50,up_scale);

[image]=fidelity(tem,im,150,up_scale);
% [image]=fidelity(im_l,aaaa,50,up_scale);

% image=ycbcr2rgb(uint8(image));
% up_scale=4;
% im_l =imresize(im_h,1/up_scale,'bicubic');
% 
% bic=imresize(im_l,up_scale,'bicubic');


%% Ω·π˚œ‘ æ
% clear
% addpath(genpath('quality')); 
% 
% load('Results/LRTMTV_DC.mat');
% % load('Results/LTTTV_DC.mat');
% % load('Results/RLTTR_DC.mat');
% [psnr4,rmse4, ergas4, sam4, uiqi4,ssim4,DD4,CC4] = quality_assessment(double(im_h), double(test_imh), 0, 0.5);
% [psnr,res_LRTMTV_DC]=csnr(im_h,test_imh,0,0);
% save('Results/res_LRTMTV_DC','res_LRTMTV_DC')

[psnr4,rmse4, ergas4, sam4, uiqi4,ssim4,DD4,CC4] = quality_assessment(double(im_h), double(image), 0, 1/up_scale);





