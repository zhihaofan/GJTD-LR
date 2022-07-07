clear;
clc;

addpath(genpath('RegularizedSC')); 
TR_IMG_PATH = 'Data';
dict_size   = 512;          % dictionary size
lambda      = 0.15;         % sparsity regularization
patch_size  = 5;            % image patch size
nSmp        = 7500;       % number of patches to sample
% nSmp        = 10500;       % number of patches to sample
upscale     = 2;            % upscaling factor

img_path=TR_IMG_PATH;
type='*.mat';
num_patch=nSmp;


[Xh, Xl] = ext_data(TR_IMG_PATH, type, patch_size, nSmp, upscale);
% [Xh, Xl] = rnd_smp_patch_s(TR_IMG_PATH, '*.bmp', patch_size, nSmp, upscale);


[Xh, Xl] = t_patch_pruning(Xh, Xl, 10); 

% save('train_data','Xh','Xl');

[Dh, Dl] = fddl_train_coupled_dict(Xh, Xl, 3,15);


% dict_path = ['Dictionary/D_' num2str(dict_size) '_' num2str(lambda) '_' num2str(patch_size) '_s' num2str(upscale) '3w.mat' ];
% dict_path=['3D_600_dim_class_5'];
save('./new_Dict_10iter_3.mat', 'Dh', 'Dl');
