function [Xh, Xl] = rnd_smp_patch_s(img_path, type, patch_size, num_patch, upscale)

img_dir = dir(fullfile(img_path, type));

Xh = [];
Xl = [];

img_num = length(img_dir);
nper_img = zeros(1, img_num);

for ii = 1:length(img_dir),
    im = imread(fullfile(img_path, img_dir(ii).name));
    nper_img(ii) = prod(size(im));
end

nper_img = floor(nper_img*num_patch/sum(nper_img));  % 需要100 000个块，每个图像提供的块的数量为nper_img的值

for ii = 1:img_num,
    patch_num = nper_img(ii);
    im = imread(fullfile(img_path, img_dir(ii).name));
    [H, L] = sample_patches(im, patch_size, patch_num, upscale);
    Xh = [Xh, H]; 
    Xl = [Xl, L];
end