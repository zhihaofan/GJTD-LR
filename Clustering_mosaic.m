clc;clear;
mosaic_dir = 'D:\data\Test_Mosaic/';
mosaic_name = dir(fullfile(mosaic_dir, '*.mat'));
total_mosaic = [];
for i=1:length(mosaic_name)
    mosaic = load(fullfile(mosaic_dir, mosaic_name(i).name));
    mosaic = mosaic.mosaic;
    vector_mosaic = reshape(mosaic, size(mosaic, 1) * size(mosaic, 2), 1);
    total_mosaic = [total_mosaic, vector_mosaic];
end
total_mosaic = double(total_mosaic);
[cls_idx,vec,cls_num]  =  Clustering(total_mosaic, 3, 15);

% function [cls_idx,vec,cls_num]  =  Clustering(x, class_num, k_iter)
% 
% 
% 
% 
% 
% 
% return;