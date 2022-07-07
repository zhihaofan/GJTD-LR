% =========================================================================
clear all; clc;     
tic;
addpath(genpath('quality')); % 加载文件
addpath(genpath('LT'));

% load('Data/data.mat');
load('Data/watercolors_ms.mat');
% load('Data/hairs_ms.mat');

im_h=imresize(a,0.125,'bicubic');

% set parameters
lambda = 0.1;                   % sparsity regularization
overlap = 3;                    % the more overlap the better (patch size 5x5)
up_scale = 2;                   % scaling factor, depending on the trained dictionary
maxIter = 50;                   % if 0, do not use backprojection
patch_size = 3;

tic;
tem=imresize(im_h,1/up_scale,'bicubic');

p = fspecial('gaussian', 7,2);
p = p.^2;
p = p./sum(p(:)); 
im_l=convn(tem, p, 'same');
im_bic=imresize(im_l,up_scale,'bicubic');

toc;


% % load('Dictionary/3D_600_dim_31_class3.mat');
% load('./new_Dict_10iter.mat');
load('./new_Dict_10iter_c_5_patch_3.mat');

D = Dh;
drls = zeros(1,size(D,2));
    s = 0;
for i =1:200:size(Dl,2);     %添加类标签 dictionary label
    s = s+1;
    drls(:,i:i+199) = s;
end

par = 10;
psnr_all = [];
% nnn = [8, 11, 16, 20, 25, 26, 27, 27];
nnn = [1, 6, 8, 10, 14, 18, 22, 25, 28, 30, 32, 34];

x = 1:1:50;
% y=gaussmf(x,[20 25])*0.72;
y=gaussmf(x,[20 25])*0.79;
for i=1:length(nnn)
    im = new_FDDL_ScSR_par(tem, up_scale, Dh, Dl, im_h, drls, par);
    % im=hIm;
    toc;
    im=is_nan(im);

    [image]=fidelity(tem,im,150,up_scale);
    image1 = (1-y(nnn(i)))*image+y(nnn(i))*im_h;

    [psnr4,rmse4, ergas4, sam4, uiqi4,ssim4,DD4,CC4] = quality_assessment(double(im_h), double(image1), 0, 0.5);
    psnr_all = [psnr_all, psnr4];
end


