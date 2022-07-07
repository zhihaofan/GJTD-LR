function [Xh, Xl] = rnd_smp_patch(img_path, type, patch_size, num_patch, upscale)

img_dir = dir(fullfile(img_path, type));

% Xh = zeros(patch_size^2,num_patch,3);
% Xl = zeros(4*patch_size^2,num_patch,3);
% Xh = zeros(patch_size^2,49969,3);      % 99966\9996
% Xl = zeros(4*patch_size^2,49969,3);

% Xh = zeros(patch_size^2,9968,3);
% Xl = zeros(4*patch_size^2,9968,3);
Xh = zeros(patch_size^2,29965,3);
Xl = zeros(4*patch_size^2,29965,3);

img_num = length(img_dir);
nper_img = zeros(1, img_num);

for ii = 1:length(img_dir)
    im = imread(fullfile(img_path, img_dir(ii).name));
    nper_img(ii) = prod(size(im));
end

[~,~,k1]=size(im);
nper_img = floor(nper_img*num_patch/sum(nper_img));  % 需要100 000个块，每个图像提供的块的数量为nper_img的值

for j=1:k1 % 光谱循环
    XH=[];
    XL=[];
    for ii = 1:img_num
        patch_num = nper_img(ii);
        im = imread(fullfile(img_path, img_dir(ii).name));
        im = im(:,:,j);
        [H, L] = sample_patches(im, patch_size, patch_num, upscale);
        XH = [XH, H]; 
        XL = [XL, L];
    end
    Xh(:,:,j)=XH;
    Xl(:,:,j)=XL;
end

% patch_path = ['Training/rnd_patches_' num2str(patch_size) '_' num2str(num_patch) '_s' num2str(upscale) '.mat'];
% save(patch_path, 'Xh', 'Xl');
