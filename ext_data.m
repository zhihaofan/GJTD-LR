function [Xh,Xl]=ext_data(img_path, type, patch_size, num_patch, upscale)

img_dir = dir(fullfile(img_path, type));
img_num = length(img_dir);
nper_img = zeros(1, img_num);

% for ii = 1:length(img_dir)
%     im = imread(fullfile(img_path, img_dir(ii).name));
%     nper_img(ii) = prod(size(im));
% end

num=length(img_dir);
num_spec=31;
a=dir(fullfile('Data/', '*.mat'));
% Xh = zeros(5^2,num_patch,num_spec);
% Xl = zeros(4*5^2,num_patch,num_spec);
Xh = zeros(patch_size^2,num_patch,num_spec);
Xl = zeros(4*patch_size^2,num_patch,num_spec);
for j=1:31 % 光谱循环
    XH=[];
    XL=[];
    for ii = 1:num
        im = load(fullfile('Data/', a(ii).name));
        temp = im.a(:,:,j);
        [H, L] = sample_patches(temp, patch_size, 1500, upscale);
        XH = [XH, H];
        XL = [XL, L];
    end
    Xh(:,:,j)=XH;
    Xl(:,:,j)=XL;
end



% num_patch=1500;
% num=1;
% num_spec=31;
% % a=dir(fullfile('Data/', '*.mat'));
% Xh = zeros(5^2,num_patch,num_spec);
% Xl = zeros(4*5^2,num_patch,num_spec);
% for j=1:31 % 光谱循环
%     XH=[];
%     XL=[];
%     for ii = 1:num
%         im = load(fullfile('Data/chart_and_stuffed_toy_ms'));
%         temp = im.a(:,:,j);
%         [H, L] = sample_patches(temp, patch_size, 1500, upscale);
%         XH = [XH, H];
%         XL = [XL, L];
%     end
%     Xh(:,:,j)=XH;
%     Xl(:,:,j)=XL;
% end









