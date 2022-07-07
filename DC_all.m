% =========================================================================
clear all; clc;     
tic;
addpath(genpath('quality')); % 加载文件
addpath(genpath('LT'));

% load('Data/data.mat');
load('Data/DC.mat');

% a=imresize(a,0.5,'bicubic');
% t=imresize(a,4,'bicubic');
% save('Data/aa','a')
% aa=aa.aa;

% a=a(8:247,8:247,:);

% a=guiyihua(a,255,0);
% im_h=aa;
im_h=DC;
% im_h=imresize(aa,0.5,'bicubic');

% set parameters
lambda = 0.1;                   % sparsity regularization
overlap = 3;                    % the more overlap the better (patch size 5x5)
up_scale = 2;                   % scaling factor, depending on the trained dictionary
maxIter = 50;                   % if 0, do not use backprojection
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
for i =1:200:size(Dl,2);     %添加类标签 dictionary label
    s = s+1;
    drls(:,i:i+199) = s;
end


im = zeros(size(im_h));
[h, w, c] = size(tem);
im_patch = 128;
step = im_patch / up_scale;
stride = 16;
for cc =1:31:191
    if cc + 31 > 191
        cc = 191 - 30;
    end
%     temp_img = DC(:, :, cc:cc+30);
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
            tem_patch = tem(x_begin + 1:x_begin+step, y_begin + 1:y_begin+step,cc:cc+30);
            im_h_patch = im(x_begin*up_scale + 1:x_begin*up_scale + im_patch, y_begin*up_scale+1:y_begin*up_scale + im_patch, cc:cc+30);
            im(x_begin*up_scale + 1:x_begin*up_scale + im_patch, y_begin*up_scale+1:y_begin*up_scale + im_patch, cc:cc+30)...
                = new_FDDL_ScSR(tem_patch, up_scale, Dh, Dl, im_h_patch, drls);
        end
    end
    
end



% im=hIm;
toc;
im=is_nan(im);


[image]=fidelity(tem,im,150,up_scale);



[psnr4,rmse4, ergas4, sam4, uiqi4,ssim4,DD4,CC4] = quality_assessment(double(im_h), double(image), 0, 0.5);

save('./Gen_data/gen_DC.mat', 'image')

