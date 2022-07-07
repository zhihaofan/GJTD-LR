function [Dh, Dl] = fddl_train_coupled_dict(Xh, Xl,class_num,k_iter)

[XX,label,cls_num] = patch_class(Xh, Xl,class_num,k_iter);
hDim=size(Xh,1);
% [D] = dictionary_training(XX, x_size, label,lambda,n_iters,cls_num);
opts.nClass        =   cls_num;
opts.wayInit       =   'random';
% opts.wayInit       =   'pca';
opts.lambda1       =   0.005;
% opts.lambda1       =   0.15;
opts.lambda2       =   0.05;
opts.nIter         =   15;
opts.show          =   true;
tr_data             =  XX;

% [Dict,Drls,coef] = FDDL(tr_data,label,opts);
[Dict,Drls,coef] = new_FDDL_SR(tr_data,label,opts);

Dh = real(Dict(1:hDim,:,:));
Dl = real(Dict(hDim+1:end, :,:));
% patch_size = sqrt(size(Dh, 1));
% dict_path = ['E:\Matlab2016a\Aworkplace\2020classify_SR\Dictionary\Dh&l_last' num2str(dict_size) '_' num2str(lambda) '_' num2str(patch_size) '_s' num2str(upscale) '.mat' ];
% save(dict_path, 'Dh', 'Dl');
